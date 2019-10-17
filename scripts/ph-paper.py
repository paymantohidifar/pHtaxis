def main():
    import funcs as mf

    # process hits obtained from blastp based on their protein domain
    mf.process_hits()

    # align TM2-proximal region of processed hits
    mf.align_hits()

    # process aligned sequences
    acc_set = mf.process_aln_seqs()

    # process taxonomy lineages of species potentially capable of pH-sensing
    mf.process_tax_lineage(acc_set)

    # remove intermediate files
    mf.clear_intermediate_files()


# ---------------------------
if __name__ == '__main__':
    main()
    print('END.')
