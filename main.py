import math as mt
from matplotlib import pyplot as plt
from abeles import Mirror


if __name__ == '__main__':
    mirror = Mirror(mt.pi / 6,  # 30 degrees
                    1.457,  # refractive index of quartz
                    2.4,  # refractive index TiO2
                    1.457)  # refractive index of quartz
    precision = 0.002
    lambda_const = 632
    number = 0  # number of double (HL) layers
    lambda_min = 400
    lambda_max = 900
    limit = 10

    # counting number of layers in mirror to achieve desired transmittance level
    mirror.set_layers()
    mirror.first_layer()
    mirror.transmittance()
    print('Dependence of mirror transmittance on its structure for TM-polarized rays:')
    print(f'{mirror.tt_tm} - AHG')
    while mirror.tt_tm > precision and number < limit:
        mirror.next_layers()
        mirror.transmittance()
        number += 1
        print(f'{mirror.tt_tm} - A(HL)^({number})HG')
    if mirror.tt_tm > precision:
        print(f'Number of double (HL) layers > {limit}, increase the limit!')
    else:
        print(f'Number of double (HL) layers = {number}')

    # calculating spectrum dots for the graph
    x = []
    y1 = []
    y2 = []
    lambda_var = lambda_min
    while lambda_var <= lambda_max:
        mirror.set_layers(lambda_const, lambda_var)
        mirror.first_layer()
        for i in range(0, number, 1):
            mirror.next_layers()
        mirror.transmittance()
        if mirror.tt_tm <= precision:
            print(f'T_TM = {mirror.tt_tm} <= {precision}, mirror is applicable for lambda = {lambda_var}')
        elif mirror.tt_tm > 1:
            print(f'T_TM(lambda = {lambda_var}) = {mirror.tt_tm} > 1, it\'s physically impossible!')
        x += [lambda_var]
        y1 += [mirror.tt_tm]
        y2 += [mirror.tt_te]
        lambda_var += 0.1

    # creating graph
    dpi = 300
    plt.figure(dpi=dpi, figsize=(1920 / dpi, 1200 / dpi))
    plt.grid(True)
    plt.xlim(lambda_min, lambda_max)
    plt.ylim(0, 1)
    plt.xlabel(r'$\lambda,\ nm$')
    plt.ylabel(r'$T$')
    plt.title("Transmittance of wavelength")
    plt.plot(x, y1, linewidth=1.0, color='Blue')
    plt.plot(x, y2, linewidth=1.0, color='Red', linestyle=':')
    plt.plot([lambda_const, lambda_const], [0, 1], linewidth=1.0, color='Green')
    plt.legend((r'$T_{TM}$', r'$T_{TE}$', r'$\lambda=632\ nm$'))
    plt.savefig('spectrum.png')
