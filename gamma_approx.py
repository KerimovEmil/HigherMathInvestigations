import math
from bernoulli_numbers import BernoulliNumber
import decimal
decimal.getcontext().prec = 100

pi = decimal.Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211"
                     "7067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954"
                     "9303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072"
                     "6024914127372458700660631558817488152092096282925409171536436789259036001133053054882046652138414"
                     "6951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891"
                     "2279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523"
                     "8467481846766940513200056812714526356082778577134275778960917363717872146844090122495343014654958"
                     "5371050792279689258923542019956112129021960864034418159813629774771309960518707211349999998372978"
                     "0499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865"
                     "8753320838142061717766914730359825349042875546873115956286388235378759375195778185778053217122680"
                     "661300192787661119590921642019893809525720106548586327")

e = decimal.Decimal("2.718281828459045235360287471352662497757247093699959574966967627724076630353547594571382178525166"
                    "42742746639193200305992181741359662904357290033429526059563073813232862794349076323382988075319525"
                    "10190115738341879307021540891499348841675092447614606680822648001684774118537423454424371075390777"
                    "44992069551702761838606261331384583000752044933826560297606737113200709328709127443747047230696977"
                    "20931014169283681902551510865746377211125238978442505695369677078544996996794686445490598793163688"
                    "92300987931277361782154249992295763514822082698951936680331825288693984964651058209392398294887933"
                    "20362509443117301238197068416140397019837679320683282376464804295311802328782509819455815301756717"
                    "36133206981125099618188159304169035159888851934580727386673858942287922849989208680582574927961048"
                    "41984443634632449684875602336248270419786232090021609902353043699418491463140934317381436405462531"
                    "52096183690888707016768396424378140592714563549061303107208510383750510115747704171898610687396965"
                    "5212671546889570350354")

gamma = decimal.Decimal('0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329'
                        '1746749514631447249807082480960504014486542836224173997644923536253500333742937337737673942792'
                        '5952582470949160087352039481656708532331517766115286211995015079847937450857057400299213547861'
                        '4669402960432542151905877553526733139925401296742051375413954911168510280798423487758720503843'
                        '1093997361372553060889331267600172479537836759271351577226102734929139407984301034177717780881'
                        '5495706610750101619166334015227893586796549725203621287922655595366962817638879272680132431010'
                        '4765059637039473949576389065729679296010090151251959509222435014093498712282479497471956469763'
                        '1850667612906381105182419744486783638086174945516989279230187739107294578155431600500218284409'
                        '6053772434203285478367015177394398700302370339518328690001558193988042707411542227819716523011'
                        '0735658339673487176504919418123000406546931429992977795693031005030863034185698032310836916400'
                        '258929708909854868257773642882539549258736295961332985747393023734388470703702844129')

d_2 = decimal.Decimal('2')
d_4 = decimal.Decimal('4')
d_5 = decimal.Decimal('5')
d_6 = decimal.Decimal('6')
d_8 = decimal.Decimal('8')
d_9 = decimal.Decimal('9')
d_11 = decimal.Decimal('11')
d_30 = decimal.Decimal('30')
d_240 = decimal.Decimal('240')
d_1620 = decimal.Decimal('1620')

lnpi = pi.ln()
ln2 = d_2.ln()


def sinh(x):
    return (e**x - e**(-x)) / 2


class Common:
    def __init__(self, x):
        self.x = x
        self.lnx = x.ln()

        self.c1 = lnpi/2 + x * (self.lnx - 1)
        self.c2 = self.c1 + (ln2 + self.lnx) / 2  # ln(c3)
        self.c3 = (2*pi*x).sqrt() * (x/e)**x  # e^(c2)


class RamanujanType:
    def __init__(self, x, c_obj):
        self.x = x
        self.common = pi.sqrt() * (x/e)**x
        self.ln_common = lnpi / 2 + x * (c_obj.lnx - 1)

    def lb(self):
        """Lower bound"""
        x = self.x
        return self.common * (8 * x ** 3 + 4 * x ** 2 + x + 1 / d_30 - 11 / (240 * x) + 5 / (240 * x ** 2))**(1/d_6)

    def ub(self):
        """Upper bound"""
        x = self.x
        return self.common * (8 * x ** 3 + 4 * x ** 2 + x + 1 / d_30 - 11 / (240 * x) + 9 / (240 * x ** 2))**(1/d_6)

    def ln_lb(self):
        """ln lower bound"""
        x = self.x
        return self.ln_common + 1/d_6 * (8 * x**3 + 4 * x**2 + x + 1/d_30 - 11/(240*x) + 5 / (240 * x**2)).ln()

    def ln_ub(self):
        """ln upper bound"""
        x = self.x
        return self.ln_common + 1/d_6 * (8 * x**3 + 4 * x**2 + x + 1/d_30 - 11/(240*x) + 9 / (240 * x**2)).ln()

    def __repr__(self):
        return 'Ramanujan type approximation, ln(pi)/2 + x(lnx-1) + 1/6 ln(8x^3 + 4x^2 + x - 11/(240x) +(5 or 9)/(240x^2))'


class SinhApproximation:
    def __init__(self, x, c_obj):
        self.x = x
        self.c_obj = c_obj

        self.ln_common = self.c_obj.c2 + (x/2) * (c_obj.lnx + sinh(1/x).ln())
        self.common = c_obj.c3 * (x*sinh(1/x))**(x/2)

    def lb(self):
        """Lower bound"""
        return self.common

    def ub(self):
        """Upper bound"""
        x = self.x
        beta = 1 / d_1620
        return self.common * (1 + beta*(x**-5))

    def ln_lb(self):
        """ln lower bound"""
        return self.ln_common

    def ln_ub(self):
        """ln upper bound"""
        x = self.x
        beta = 1 / d_1620
        return self.ln_common + (1 + beta*(x**-5)).ln()

    def __repr__(self):
        return 'Sinh approximation, (ln(2) + ln(pi)+ ln(x))/2 + x(lnx-1) + x/2 *(ln(x) + ln(sinh(1/x)) OR -5ln(x)-ln(1620)'


class RobbinsApproximation:
    def __init__(self, x, c_obj):
        self.x = x
        self.c_obj = c_obj

    def lb(self):
        """Lower bound"""
        x = self.x
        return self.c_obj.c3 * e**(1/(12*x+1))

    def ub(self):
        """Upper bound"""
        x = self.x
        return self.c_obj.c3 * e ** (1 / (12 * x))

    def ln_lb(self):
        """ln lower bound"""
        return self.c_obj.c2 + 1/(12*self.x + 1)

    def ln_ub(self):
        """ln upper bound"""
        return self.c_obj.c2 + 1/(12*self.x)

    def __repr__(self):
        return 'Robbins approximation, (ln(2) + ln(pi)+ ln(x))/2 + x(lnx-1) + 1/(12x + (0 OR 1))'


class XinLi_ChaoPingChenApproximation:
    """"https://www.emis.de/journals/JIPAM/images/264_06_JIPAM/264_06.pdf"""
    def __init__(self, x, c_obj):
        self.x = x
        self.c_obj = c_obj

        self.common = ((x+1)/e)**(x+1) * e
        # self.ln_common = x*self.c_obj.lnx - x + 1
        self.ln_x1 = (x+1).ln()
        self.ln_common = (x+1)*self.ln_x1 - x

    def lb(self):
        """Lower bound"""
        return self.common * (self.x + 1)**(-gamma)

    def ub(self):
        """Upper bound"""
        return self.common * (self.x + 1)**(-1/d_2)

    def ln_lb(self):
        """ln lower bound"""
        return self.ln_common - gamma*self.ln_x1

    def ln_ub(self):
        """ln upper bound"""
        return self.ln_common - self.ln_x1/d_2

    def __repr__(self):
        return 'Xin Li and Chao-Ping Chen approximation, x^(x-gamma) < Gamma(x)*e^(x-1) < x^(x-1/2)'


class BernoulliType:
    def __init__(self, x, c_obj, bern_num=7):
        self.x = x
        self.common = c_obj.c3
        self.ln_common = c_obj.c2

        bernoulli_obj = BernoulliNumber()
        self.ls_b = [decimal.Decimal(float(
            bernoulli_obj.get(2*i)/(2*i*(2*i-1))))/(x**(2*i - 1))
                     for i in range(1, bern_num)]

        # since the bernoulli series is alternating
        # if odd bernoulli number given, then this is an upper bound
        # if even number then this is a lower bound

        if bern_num % 2 == 1:  # even number, then this is a upper bound
            self.lower_b = sum(self.ls_b)
            self.upper_b = sum(self.ls_b[:-1])
        else:
            self.lower_b = sum(self.ls_b[:-1])
            self.upper_b = sum(self.ls_b)

    def lb(self):
        """Lower bound"""
        return self.common * e ** self.lower_b

    def ub(self):
        """Upper bound"""
        return self.common * e ** self.upper_b

    def ln_lb(self):
        """ln lower bound"""
        return self.ln_common + self.lower_b

    def ln_ub(self):
        """ln upper bound"""
        return self.ln_common + self.upper_b

    def __repr__(self):
        return 'Bernoulli type approximation, ln(n)(n + 1/2) - n -ln(2pi) + (1/12n - 1/360n^3 + 1/1260n^5 + O(1/n^7)'


def f_print(n, lb, a, ub):
    print("n:{:.5}, lower:{}, actual:{}, higher:{}".format(n, lb, a, ub))
    if lb >= a or ub <= a:
        raise ValueError(lb, a, ub)
    else:
        print("lb_diff:{:.2e}, ub_diff:{:.2e}".format(a - lb, ub - a))
        print("lb_diff_pct:{:.2e}, ub_diff_pct:{:.2e}".format(1 - lb / a, ub / a - 1))
    print('----------------------------------------------------------------------------')


def approximation_gamma(max_n, ls_approx):
    # Log gamma
    for i in range(2, max_n+1):
        # context
        x = decimal.Decimal(i)
        common_obj = Common(x)

        # Actual
        actual = decimal.Decimal(math.factorial(x))

        print('----------------------------------------------------------------------------')
        print("x = {}, actual = {}".format(x, actual))
        print('----------------------------------------------------------------------------')

        # Approximations
        for approx in ls_approx:
            approx_obj = approx(x=x, c_obj=common_obj)
            print(approx_obj)
            f_print(x, approx_obj.lb(), actual, approx_obj.ub())


def approximation_ln_gamma(max_n, ls_approx):
    # Log gamma
    for i in range(2, max_n+1):
        # context
        x = decimal.Decimal(i)
        common_obj = Common(x)

        # Actual
        actual = decimal.Decimal(math.factorial(i)).ln()

        print('----------------------------------------------------------------------------')
        print("x = {}, actual = {}".format(x, actual))
        print('----------------------------------------------------------------------------')

        # Approximations
        for approx in ls_approx:
            approx_obj = approx(x=x, c_obj=common_obj)
            print(approx_obj)
            f_print(x, approx_obj.ln_lb(), actual, approx_obj.ln_ub())


if __name__ == '__main__':
    ls_approximations = [RamanujanType, SinhApproximation, RobbinsApproximation, XinLi_ChaoPingChenApproximation,
                         BernoulliType]
    approximation_ln_gamma(max_n=20, ls_approx=ls_approximations)
    # approximation_gamma(max_n=20, ls_approx=ls_approximations)
