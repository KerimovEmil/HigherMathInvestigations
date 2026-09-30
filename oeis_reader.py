import requests
import json


def load_oeis_sequence_table(sid, max_n=None):
    r"""
    Gets the table of terms of the sequence #sid (should regex match
    r'A(\d+) (.*)') from its remote b-file, e.g. the b-file for the
    sequence A001221 is located at
        https://oeis.org/A003415/b003415.txt
    More information can be found here:
        https://oeis.org/wiki/B-files
    """
    sid = sid.upper()
    num_part = sid[1:] if sid.startswith('A') else sid
    res = requests.get(f'https://oeis.org/A{num_part}/b{num_part}.txt')
    res.raise_for_status()

    table = []
    for line in res.text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                table.append(int(parts[1]))
            except ValueError:
                continue

    if max_n is not None:
        table = table[:max_n]

    return table


def get_oeis_sequence_meta(sid, key='name'):
    """
    Returns the sequence metadata dictionary - more information can be
    found here:
        https://oeis.org/wiki/JSON_Format,_Compressed_Files
    key: 'name', 'comment'
    """
    sid = sid.upper()
    res = requests.get(f'https://oeis.org/search?q={sid}&fmt=json')
    res.raise_for_status()
    meta_data = res.json()

    if meta_data and meta_data.get('results'):
        return meta_data['results'][0].get(key)
    return None


if __name__ == '__main__':
    print(load_oeis_sequence_table("A003415", 10))
    print(get_oeis_sequence_meta("A003415"))
