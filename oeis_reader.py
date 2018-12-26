import requests


def load_oeis_sequence_table(sid, max_n=None):
    """
        Gets the table of terms of the sequence #sid (should regex match
        'A(\d+) (.*)') from its remote b-file, e.g. the b-file for the
        sequence A001221 is located at
            https://oeis.org/A003415/b003415.txt
        More information can be found here:
            http://oeis.org/wiki/B-files
    """
    res = requests.get('http://oeis.org/{}/b{}.txt'.format(sid, sid[1:]))
    table = [int(s.split(' ')[1]) for s in res.text.split('\n') if ' ' in s]
    if max_n is not None:
        table = table[:max_n]

    return table


if __name__ == '__main__':
    print(load_oeis_sequence_table("A003415", 10))
