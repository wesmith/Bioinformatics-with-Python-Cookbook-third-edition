# ws_utils.py
# WESmith 11/10/22


def attrs(obj, skip=True, token='__'):
    '''
    Convenience function to print an objects attributes.
    '''
    attr = ['OBJECT TYPE: {}'.format(type(obj))]
    for k in dir(obj):
        if skip and k.__contains__(token): continue
        attr.append(k)
    return attr


def phred_to_percent_accurate(phred):
    return 100 * (1. - 10**(-phred/10))