"""An interesting way to divide"""


def new_mod(a, m):
    """
    Returns a mod m.
    Works well for m=0,1,2,3,4,5,8,9,10,11
    Args:
        a: <int>
        m: <int>
    Returns: a mod m
    """
    str_a = str(a)
    if len(str_a) > 2:

        if m == 0 or m == 1 or a == m:
            return 0

        if m == 2:
            last = int(str_a[-1:])
            return new_mod(last, m)

        if m == 3 or m == 9:
            sum_of_digits = sum([int(d) for d in str_a])
            return new_mod(sum_of_digits, m)

        if m == 4:
            last = int(str_a[-1])
            second_last = int(str_a[-2:-1])
            answer = 2 * second_last + last
            return new_mod(answer, m)

        if m == 5:
            last = int(str_a[-1])
            return new_mod(last, m)

        if m == 7:
            last = int(str_a[-1:])
            first = int(str_a[:-1])
            answer = new_mod(first - 2 * last, m)
            if answer == 0:
                return 0
            else:
                return a % m

        if m == 8:
            last = int(str_a[-1:])
            second_last = int(str_a[-2:-1])
            third_last = int(str_a[-3:-2])
            answer = 4 * third_last + 2 * second_last + last
            return new_mod(answer, m)

        if m == 10:
            last = int(str_a[-1:])
            return last

        if m == 11:
            new_a = 0
            for i, digit in enumerate(str_a):
                if not i % 2:
                    new_a += int(digit)
                else:
                    new_a -= int(digit)
            return new_mod(new_a, m)

        if m == 13:
            last = int(str_a[-1:])
            first = int(str_a[:-1])
            answer = new_mod(first - 9 * last, m)
            if answer == 0:
                return 0
            else:
                return a % m

        return a % m

    else:

        return a % m


if __name__ == '__main__':
    assert new_mod(400, 3) == 1

    test_num = 442332523255252342423453323536236363463246346111111111422
    assert new_mod(test_num, 2) == 0
    assert new_mod(test_num, 3) == 2
    assert new_mod(test_num, 4) == 2
    assert new_mod(test_num, 5) == 2
    assert new_mod(test_num, 6) == 2
    assert new_mod(test_num, 7) == 1
    assert new_mod(test_num, 8) == 6
    assert new_mod(test_num, 9) == 8
    assert new_mod(test_num, 10) == 2
    assert new_mod(test_num, 11) == 0
    assert new_mod(test_num, 13) == 12
