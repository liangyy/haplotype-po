def read_table(fname):
    if os.path.splitext(fname) == 'gz':
        kw = {'compression' : 'gzip'}
    else:
        kw = {}
    return pd.read_csv(fname, sep='\t', **kw)