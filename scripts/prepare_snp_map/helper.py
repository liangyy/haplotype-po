import pandas as pd
def make_id(x):
    return '{chr_}:{pos}:{a1}:{a2}'.format(chr_=x[0], pos=x[1], a1=x[2], a2=x[3])
def join_with_varid(candidate_id1, candidate_id2, id_pool):
    df = pd.DataFrame({
        'id1': candidate_id1,
        'id2': candidate_id2
    })
    pool = pd.Series(id_pool)
    df['id1_check'] = df['id1'].isin(pool)
    df['id2_check'] = df['id2'].isin(pool)
    df['check_combine'] = df[['id1_check', 'id2_check']].apply(lambda x: x.sum(), axis=1)
    df['assigned_id'] = df[['id1', 'id2', 'id1_check', 'id2_check', 'check_combine']].apply(lambda x: _assign_id(x), axis=1)
    return df
def _assign_id(x):
    if x.check_combine == 2:
        return 'ambiguious'
    if x.check_combine == 0:
        return 'not_shown'
    if x.id1_check is True:
        return x.id1
    if x.id2_check is True:
        return x.id2