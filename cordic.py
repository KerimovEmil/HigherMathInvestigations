# https://en.wikipedia.org/wiki/CORDIC

# v_0 = [ 1 ]
#       [ 0 ]

# v_i = R_i * v_{i-1}

# R_i = [ cos(gamma_i) -sin(gamma_i) ]
#       [ sin(gamma_i) cos(gamma_i)  ]

# R_i = 1/(sqrt(1 + tan^2(gamma_i))) * [ 1            -tan(gamma_i) ]
#                                      [ tan(gamma_i)       1       ]

# pick gamma_i such that tan(gamma_i) = +/- 2^(-i)

if __name__ == '__main__':
    n = 30
    # K
    k = 1
    # see how fast this product converges
    for i in range(n):
        print(k)
        k *= 1 / (1 + 2**(-2*i))**0.5
    print(k)
