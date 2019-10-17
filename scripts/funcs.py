import os
import matplotlib.pyplot as plt
from re import split


def read_fasta_file(filename):
    """ Reads in file in fasta format and returns sequences of tuples containing
    header name and the corresponding sequence"""
    with open(filename, 'r') as file:
        return tuple((name, seq.replace('\n', ''))
                     for name, ignore, seq in [entry.partition('\n')
                                               for entry in file.read().split('>')[1:]])


def blastp_parser(filename):
    """ Reads in blastp output file (-outfmt 6 & 7 ) and returns content of table as sequences of
    keys and values, where keys holds the attributes in the table"""
    with open(filename, 'r') as file:
        table = file.read().split('#')
        attribs = table[4].replace(' Fields: ', '').rstrip('\n').split(', ')
        return tuple({attribs[n]: value for n, value in enumerate(hit.split('\t'))}
                     for hit in table[5].split('\n')[1:-1])


def hmmscan_parser(filename):
    """ Reads in hmmscan output file (tblout format) and returns tables attributes as
     sequence of keys and values. Note that only the best domain for each query is kept."""
    with open(filename, 'r') as file:
        domains = tuple()
        for line in file.read().split('\n'):
            if line.startswith('#') or line == '':
                continue
            line = split('\s+', line, 22)
            domain = dict()
            domain.setdefault('target name', line[0])
            domain.setdefault('accession', line[1])
            domain.setdefault('query name', line[2])
            domain.setdefault('domain evalue', line[7])
            domain.setdefault('description', ' '.join([line[n] for n in range(18, len(line))]))
            domains += (domain,)
        domains = collapse_duplicates(domains, 'query name')
    return domains


def alignment_file_parser(filename):
    """ Reads alignment file (in fasta format) and returns sequence of header name
      and the corresponding sequence. Here is alignment file format:
    >accession number|species name
    sequence
    """
    aln_hits = tuple()
    for entry in read_fasta_file(filename):
        hit = dict()
        hit['subject acc.'] = entry[0].split('|')[0]
        hit['species name'] = entry[0].split('|')[1]
        hit['subject seq'] = entry[1]
        aln_hits += (hit,)
    return aln_hits


def get_species_taxon(hits):
    """ Takes in hits and returns hits which belong to bacterial species taxon."""
    with open(dbdir + '2.species.taxids', 'r') as file:
        taxids = {line.rstrip('\n') for line in file.readlines()}
        return tuple(hit for hit in hits if hit.get('subject tax ids') in taxids)


def extract_seqs(hits, start, end=None):
    """ Takes in hits and the desired range (start, end) and update sequence field
    for all hits"""
    for n in range(0, len(hits)):
        hits[n]['subject seq'] = hits[n]['subject seq'][start:end]


def collapse_duplicates(hits, key):
    """ Collapses hits with 100% similar sequences into one hit based on hits order"""
    duplicates = set()
    return tuple(hit for hit in hits
                 if not (hit.get(key) in duplicates or duplicates.add(hit.get(key))))


def get_redundancy_level(seq1, seq2):
    """ Calculates similarity levels of two sequences"""
    return list(map(lambda x: x[0] == x[1], list(zip(seq1, seq2)))).count(True) / len(seq1)


def collapse_redundant_hits(hits, threshold=0.9):
    """ Collapses hits with similarity level set by threshold argument. Default
    threshold is 90%"""
    grouped_hits = []
    while len(hits) > 1:
        prev_hit = hits[0]
        subgroup = [prev_hit]
        ungroup = []
        for hit in hits[1:]:
            if get_redundancy_level(prev_hit.get('subject seq'), hit.get('subject seq')) > threshold:
                subgroup.append(hit)
            else:
                ungroup.append(hit)
        grouped_hits.append(subgroup)
        hits = ungroup
    return tuple(subgroup[0] for subgroup in grouped_hits)


def collect_loci(hits, upper_locus_range, lower_locus_range):
    """ Takes in hits and updates upper locus and lower locus sequence field using loci ranges
     provided as arguments """
    for n in range(0, len(hits)):
        hits[n]['upper locus'] = hits[n]['subject seq'][upper_locus_range[0]:upper_locus_range[1]]
        hits[n]['lower locus'] = hits[n]['subject seq'][lower_locus_range[0]:lower_locus_range[1]]
        hits[n]['analyzed'] = False


def locus_type(hits, upper_locus_type, lower_locus_type):
    """ Takes in hits and returns hits whose upper and lowe loci match the provided loci types
    as arguements"""
    aa_set = {'charged': 'RHKDE', 'polar': 'STYNQ', 'nonpolar': 'GAVCPLIMWF', 'not_aligned': '-'}
    selected = tuple()
    for n in range(0, len(hits)):
        if (len(set(hits[n]['upper locus']) & set(aa_set[upper_locus_type])) > 0 and
                len(set(hits[n]['lower locus']) & set(aa_set[lower_locus_type])) > 0 and not hits[n]['analyzed']):
            selected += (hits[n],)
            hits[n]['analyzed'] = True
    return selected


def pie_plot(data, save_path, show_fig=False):
    """ generates pie plot"""

    def set_wedge_numbers(perc):
        return '{:.1f}%'.format(perc)

    plt.pie(data, explode=(0.02, 0.02, 0.02), labels=['charged-charged', 'charged-polar', 'other'],
            colors=['green', 'orange', 'gray'], autopct=lambda perc: set_wedge_numbers(perc), center=(0, 0),
            textprops=dict(color='black', fontsize=10))
    plt.savefig(save_path, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None,
                format='eps', transparent=False, bbox_inches=None, pad_inches=0.1, metadata=None)
    if show_fig:
        plt.show()


def process_domains():
    """ Runs hmmscan to obtain the best domain for each hit"""
    # run hmmcan to identify the best domain matches for each hit
    print('Running hmmscan to search each hit against Pfam 32.0 hmm profile database to determine the best domain '
          'matched\n(May take a few minutes ...)')
    os.system('hmmscan --cut_ga --cpu 4 --tblout {} {} {} > {}'
              .format(resultdir + 'domains.tab', pfamdb, resultdir + '_processed_hits.fa', logsdir + 'hmmscan.log'))

    domains = hmmscan_parser(resultdir + 'domains.tab')
    print('Total {} of unique domains were parsed for each hit.'.format(len(domains)))

    print('Processing domains ...')
    domain_types = {domain['target name'] for domain in domains}
    for domain_type in domain_types:
        evalues = tuple(float(domain['domain evalue']) for domain in domains if domain['target name'] == domain_type)
        print('{:.2f}% of the hits have \'{}\' domain with {} < E-values < {}.'
              .format(len(evalues) / len(domains) * 100, domain_type, min(evalues), max(evalues)))

    cutoff_evalue = 1e-10
    domains = tuple(domain['query name'] for domain in domains if float(domain['domain evalue']) < cutoff_evalue)
    print('Total {} of hits were remained after removing hits with E-value > {}'.format(len(domains), cutoff_evalue))
    return domains


def process_hits():
    """ Processes blastp output """
    # parse hit table obtained form blastp search
    print('Parsing hits obtained from blastp search ...')
    hits = blastp_parser(datadir + '2.species.hits.tab')
    print('Total {} of hits were parsed from blastp search output.'.format(len(hits)))

    print('Processing the hits ...')
    # restrict the hits to species taxonomy level
    hits = get_species_taxon(hits)
    print('Total {} of hits remained after restricting hits to only species taxonomy level.'.format(len(hits)))

    # remove gaps from sequences
    for n in range(0, len(hits)):
        hits[n]['subject seq'] = hits[n]['subject seq'].replace('-', '')

    # remove hits, which did not align to membrane-proximal module of McpA (163-239)
    hits = tuple(hit for hit in hits
                 if (int(hit['q. end']) > 239 and int(hit['q. start']) < 163))
    print('Total {} of hits remained after removing the hits that did not align to membrane-proximal module.'
          .format(len(hits)))

    # collect the hits that have 80% coverage of the query sequence (247 aa).
    hits = tuple(hit for hit in hits if len(hit['subject seq']) > 0.8*247)
    print('Total {} of hits remained after removing the hits that had less than 80% coverage [excluding gaps].'
          .format(len(hits)))

    # remove duplicate sequences
    hits = collapse_duplicates(hits, 'subject seq')
    print('Total {} of hits remained after removing duplicate sequences.'.format(len(hits)))

    # write processed hits in a fasta file used for domain identification
    with open(resultdir + '_processed_hits.fa', 'w') as file:
        file.writelines('>' + hit['subject acc.'] + '\n' + hit['subject seq'] + '\n' for hit in hits)
    print(100*'-')

    # identify domains for each hit
    domains = process_domains()

    # Collect hits that satisfied allowed threshold for domain
    hits = tuple(hit for hit in hits if hit['subject acc.'] in domains)

    # extract aa sequence of pH sensing region (membrane-proximal module) of collected hits
    extract_seqs(hits, start=-90)
    print('Sequences of TM2-proximal region ({} aa) were extracted.'.format(90))

    # remove hits that have same aa sequence to reduce bias
    hits = collapse_duplicates(hits, 'subject seq')
    print('Total {} of hits remained after removing duplicate sequences.'.format(len(hits)))

    # write processed hits in a fasta file
    with open(resultdir + '_seqs-for-align.fa', 'w') as file:
        file.writelines('>' + '|'.join([hit['subject acc.'],
                                        hit['subject title'].rstrip(']').partition('[')[2].strip('[').replace(' ', '_')])
                        + '\n' + hit['subject seq'] + '\n' for hit in hits)
    print(100*'-')


def align_hits():
    """ Runs multiple-sequence alignments on processed hits"""
    print('Running clustal omega for fast initial MSA ...')
    os.system('clustalo -i {} -o {} --force > /dev/null'
              .format(resultdir + '_seqs-for-align.fa', resultdir + '_aligned_seqs.fa'))

    print('Running trimal to remove badly aligned sequences ...')
    os.system('trimal -in {} -out {} -resoverlap 0.80 -seqoverlap 90 > /dev/null'
              .format(resultdir + '_aligned_seqs.fa', resultdir + '_aligned_seqs_trimmed.fa'))

    aln_seqs = alignment_file_parser(resultdir + '_aligned_seqs_trimmed.fa')
    print('Total {} of hits remained after removing badly aligned hits.'.format(len(aln_seqs)))

    print('Trimming alignments and collapsing highly similar hits with 95% cutoff ...')
    left_cut_index = max([aln_seq['subject seq'].find('TKKVN') for aln_seq in aln_seqs])
    right_cut_index = max([aln_seq['subject seq'].find('AAQP') + 4 for aln_seq in aln_seqs])
    extract_seqs(aln_seqs, left_cut_index, right_cut_index)
    aln_seqs = collapse_redundant_hits(aln_seqs, 0.95)
    with open(resultdir + '_aligned_seqs_trimmed.fa', 'w') as file:
        file.writelines('>' + '|'.join([hit['subject acc.'], hit['species name'].replace(' ', '_')])
                        + '\n' + hit['subject seq'].replace('-', '') + '\n' for hit in aln_seqs)
    print('Total {} of hits remained after collapsing highly similar sequences with {}% similarity.\n'
          .format(len(aln_seqs), 0.95))

    print('Running Muscle for accurate MSA ...')
    os.system('muscle -in {} -out {} -log {} > /dev/null'
              .format(resultdir + '_aligned_seqs_trimmed.fa', resultdir + 'aligned_seqs.fa', logsdir + 'muscle.logs'))


def process_aln_seqs():
    """ Processes alignment sequences"""
    aln_seqs = alignment_file_parser(resultdir + 'aligned_seqs.fa')

    # collect pH residues from each alignment and sort alignments based on species names
    upper_ph_res = max([aln_seq['subject seq'].find('TQGYAFI') for aln_seq in aln_seqs])
    lower_ph_res = max([aln_seq['subject seq'].find('HEAAQP') for aln_seq in aln_seqs])
    collect_loci(aln_seqs, (upper_ph_res, upper_ph_res + 2), (lower_ph_res, lower_ph_res + 2))
    for n in range(0, len(aln_seqs)):
        aln_seqs[n].setdefault('analyzed', False)
    cc_aln_seqs = locus_type(aln_seqs, 'charged', 'charged')
    cp_aln_seqs = locus_type(aln_seqs, 'charged', 'polar')
    pc_aln_seqs = locus_type(aln_seqs, 'polar', 'charged')
    ph_aln_seqs = cc_aln_seqs + cp_aln_seqs + pc_aln_seqs
    ph_aln_seqs = sorted(ph_aln_seqs, key=lambda x: x.get('upper locus') + x.get('lower locus'), reverse=True)
    ph_aln_seqs = sorted(ph_aln_seqs, key=lambda x: x.get('species name'), reverse=False)
    with open(resultdir + 'ph_aligned_seqs.fa', 'w') as file:
        file.writelines('>' + '|'.join([hit['subject acc.'], hit['species name']]) +
                        '\n' + hit['subject seq'] + '\n' for hit in ph_aln_seqs)
    print('Hits that are potentially capable of pH-sensing are written to <{}>.'
          .format(resultdir + 'ph_aligned_seqs.fa'))

    fractions = (len(cc_aln_seqs)/len(aln_seqs), (len(cp_aln_seqs) + len(pc_aln_seqs)) / len(aln_seqs),
                 1 - len(ph_aln_seqs) / len(aln_seqs))
    pie_plot(fractions, figuredir + 'pie.eps', True)
    print('Pie plot representing chemistry of pH amino acid residues is saved to <{}>.'.format(figuredir + 'pie.sps'))

    # collect hits which have highly conserved pH-sensing residues wrt to Bs 168 pH chemoreceptors
    with open(resultdir + 'Bs168_conserved_hits.fa', 'w') as file:
        file.writelines('>' + '|'.join([hit['subject acc.'], hit['species name']]) + '\n' +
                        hit['subject seq'].replace('-', '') + '\n' for hit in aln_seqs
                        if (hit['upper locus'] in {'KE', 'TQ'} and hit['lower locus'] in {'HE', 'HD', 'QD', 'KD'}))

    # collect upper and lower pH-sensing loci
    collect_loci(ph_aln_seqs, (upper_ph_res - 4, upper_ph_res + 6), (lower_ph_res - 4, lower_ph_res + 6))
    with open(resultdir + 'upper_loci_logo.fa', 'w') as file1, open(resultdir + 'lower_loci_logo.fa', 'w') as file2:
        file1.writelines('>' + hit['subject acc.'] + '\n' + hit['upper locus'] + '\n' for hit in ph_aln_seqs)
        file2.writelines('>' + hit['subject acc.'] + '\n' + hit['lower locus'] + '\n' for hit in ph_aln_seqs)
    print('pH loci sequences are written to <{}> and <{}>.'
          .format(resultdir + 'lower_loci_logo.fa', resultdir + 'upper_loci_logo.fa'))

    print(100*'-')
    return {hit['subject acc.'] for hit in ph_aln_seqs}


def process_tax_lineage(acc_set):
    """ Processes taxonomy lineage of identified bacterial species"""
    # extract taxonomy lineage of all identified species
    hits = blastp_parser(datadir + '2.species.hits.tab')
    acc_taxids = {(hit['subject acc.'], hit['subject tax ids']) for hit in hits}
    ph_taxids = {elt[1] for elt in acc_taxids if elt[0] in acc_set}

    print('Running taxonkit to extract taxonomy lineages for species potentially capable of pH-taxis...')
    with open(resultdir + 'ph_aligned_seqs.taxids', 'w') as file:
        file.writelines(taxid + '\n' for taxid in ph_taxids)
    os.system('taxonkit lineage {} -o {}'.format(resultdir + 'ph_aligned_seqs.taxids', resultdir + '_lineage.txt'))
    reformat = '{p};{c};{o};{f};{g};{s}'
    os.system('taxonkit reformat {} -f \'{}\' | tee {} > /dev/null'
              .format(resultdir + '_lineage.txt', reformat, resultdir + '_ph_aligned_seqs.lineage'))

    tax_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
    with open(resultdir + '_ph_aligned_seqs.lineage', 'r') as file:
        tax_entries = tuple({tax_levels[n]: value for n, value in enumerate(entry.split(';'))} for entry in
                            tuple(line.rstrip('\n').split('\t')[-1] for line in file.readlines()))
        tax_entries = sorted(tax_entries, key=lambda x: x.get('species'))

    with open(resultdir + 'ph_aligned_seqs.lineage.tsv', 'w') as file:
        file.write('\t'.join(tax_levels) + '\n')
        file.writelines('\t'.join(list(entry.values())) + '\n' for entry in tax_entries)
    print('Taxonomy lineages are written to <{}>'.format(resultdir + 'ph_aligned_seqs.lineage.tsv'))

    for level in tax_levels:
        print('Total number of unique {}: {}'.format(level, len({entry.get(level) for entry in tax_entries})))
    print(100*'-')


def clear_intermediate_files():
    """ Remove all intermediate files generated during analysis for better organization"""
    os.system('rm -f {}_*'.format(resultdir))


def main():
    """ DO NOTHING """
    pass


# ---------------------------
if __name__ == '__main__':
    main()
else:
    datadir = '../data/'
    resultdir = '../results/'
    figuredir = '../results/figures/'
    logsdir = '../logs/'
    dbdir = '../db/'
    pfamdb = '$HOME/hmmer-3.2.1/db/Pfam-A.hmm'
