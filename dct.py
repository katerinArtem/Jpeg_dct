

# import cv2
# import numpy as np
from math import cos,pi

#RGB → YCbCr → 4:2:0 subsampling → DCT → квантование → обратно.

# def demo_image():
#     """
#     Читает '0.jpg' и '0.png', применяет YCbCr + 4:2:0, сохраняет
#     'output_image.jpg'.
#     jim = cv2.imread("0.jpg")
#     im  = cv2.imread("0.png")

#     im  = cv2.cvtColor(im, cv2.COLOR_BGRA2RGB)
#     jim = cv2.cvtColor(im, cv2.COLOR_BGR2YCrCb)   # см. замечание выше

#     # ── sanity-check цвета на пикселе [0,0] ──
#     # NB: jim хранится в порядке YCrCb (OpenCV), а YCbCrtoRGB ждёт YCbCr,
#     #     поэтому 'CRGB' ниже фактически считается с перепутанными Cr/Cb.
#     #     В реконструкции пикселя (ниже) этот свап компенсируется.
#     print("conv_check")
#     print("ORGB",   list(im[0, 0]))
#     print("CRGB",   YCbCrtoRGB(list(jim[0, 0])))
#     print("OYCbCr", list(jim[0, 0]))
#     print("CYCbCr", RGBtoYCbCr(list(im[0, 0])))

#     ty, tx = 0, 150
#     ojim, oim = jim, im
#     jim = jim.tolist()

#     # subsampl_4_2_0 возвращает [Y, Cb_m, Cr_m], где по факту Cb_m хранит
#     # OpenCV-канал Cr, а Cr_m — Cb. Так что ниже YCbCr[2] = настоящий Cb,
#     # YCbCr[1] = настоящий Cr → реконструкция даёт корректный [Y, Cb, Cr].
#     YCbCr = subsampl_4_2_0(jim)

#     jim = [[
#         [
#             YCbCr[0][i][k],             # Y  (как есть)
#             YCbCr[2][i // 2][k // 2],   # Cb (из "Cr_m")
#             YCbCr[1][i // 2][k // 2],   # Cr (из "Cb_m")
#         ]
#         for k in range(len(YCbCr[0][0]))
#     ] for i in range(len(YCbCr[0]))]

#     y1, y2 = 0, 511
#     x1, x2 = 0, 511
#     for i in range(y1, y2):
#         for k in range(x1, x2):
#             im[i, k] = np.array(YCbCrtoRGB(jim[i][k]))

#     print("YCbCr after subs ", jim[0][7])
#     print("RGB   after subs ", YCbCrtoRGB(jim[ty][tx]))
#     print("origin YCbCr    ", ojim[ty, tx])
#     print("origin RGB      ", oim[ty, tx])

#     img = np.array(im[y1:y2, x1:x2], np.uint8)
#     cv2.imwrite('output_image.jpg', img)


def cust(N, i, k):
    """Базисная функция DCT-II: cos(π/N · (i + 0.5) · k)."""
    return cos((pi / N) * (i + 0.5) * k)

# ═══════════════════════════════════════════════════════════════════════════
# 1D DCT (прямое/обратное) + квантование
# ═══════════════════════════════════════════════════════════════════════════

def DCT(k, g_s):
    """
    Прямое 1D DCT-II:
        G[k] = Σ_{i=0}^{N-1} g_s[i] · cos(π/N · (i + 0.5) · k).
    k=0 — DC (постоянная составляющая ≈ среднее), k>0 — AC (частоты).
    """
    return sum(g_s[i] * cust(len(g_s), i, k) for i in range(len(g_s)))

def rDCT(u, G_s):
    """
    Обратное 1D DCT-III:
        s[u] = (2/N)·(G[0]/2 + Σ_{k=1}^{N-1} G[k]·cos(π/N·(u+0.5)·k)).
    Деление на 2 у DC — чтобы DCT и IDCT были взаимно-обратными.
    """
    N = len(G_s)
    return (2 / N) * (
        G_s[0] * (1 / 2)
        + sum(G_s[i] * cos((pi / N) * (u + 0.5) * i) for i in range(1, N))
    )

def QDCT(k, G_s, q_s):
    """Квантование: T[k] = round(G[k] / Q[k]). Тут возникают потери."""
    return round(G_s[k] / q_s[k])

def rQDCT(k, t_s, q_s):
    """Обратное квантование: G'[k] = T[k] · Q[k]."""
    return t_s[k] * q_s[k]

def d1ex(o_s, q_s):
    """Демо 1D-цикла: DCT → квант → де-квант → IDCT."""
    print("S:", o_s)
    G_s = [DCT(m, o_s) for m in range(len(o_s))]
    print("D:", G_s)
    print("Q:", q_s)
    t_s = [QDCT(k, G_s, q_s) for k in range(len(o_s))]
    print("T:", t_s)
    rt_s = [rQDCT(k, t_s, q_s) for k in range(len(o_s))]
    print("D':", rt_s)
    rtG_s = [rDCT(u, rt_s) for u in range(len(o_s))]
    print("S'", list(map(round, rtG_s)))
    rG_s = [rDCT(u, G_s) for u in range(len(o_s))]
    print("S",  list(map(round, rG_s)))


# ═══════════════════════════════════════════════════════════════════════════
# 2D DCT (сепарабельно) + квантование матрицей
# ═══════════════════════════════════════════════════════════════════════════

def MDCT(o_s, k1, k2):
    """
    2D DCT, коэффициент (k1, k2):
        G[k1,k2] = Σ_i1 Σ_i2 o_s[i1][i2] · cos_y(i1,k1) · cos_x(i2,k2).
    """
    return sum(DCT(k2, o_s[i1]) * cust(len(o_s), i1, k1)
               for i1 in range(len(o_s)))


def rMDCT(G_s, k1, k2):
    """Обратное 2D DCT (сепарабельно: IDCT по X, потом по Y)."""
    N = len(G_s)
    return (2 / N) * (
        rDCT(k2, G_s[0]) * (1 / 2)
        + sum(rDCT(k2, G_s[i1]) * cos((pi / N) * (k1 + 0.5) * i1)
              for i1 in range(1, N))
    )


def qMDCT(G_s, Q_s):
    """Поэлементное квантование 2D-матрицы коэффициентов."""
    return [[round(g / q) for g, q in zip(G_row, Q_row)]
            for G_row, Q_row in zip(G_s, Q_s)]


def rqMDCT(T_s, Q_s):
    """Поэлементное обратное квантование."""
    return [[t * q for t, q in zip(T_row, Q_row)]
            for T_row, Q_row in zip(T_s, Q_s)]

def d2ex(o_s, q_s):
    """Демо 2D-цикла: DCT → квант → де-квант → IDCT + MSE."""
    print("orig")
    for row in o_s: print(row)

    print("mdct")
    mdct = [[MDCT(o_s, k1, k2) for k2 in range(len(o_s[k1]))]
            for k1 in range(len(o_s))]
    for row in mdct: print(row)

    print("rmdct  (без квантования — проверка обратимости)")
    rmdct = [[rMDCT(mdct, k1, k2) for k2 in range(len(mdct[k1]))]
             for k1 in range(len(mdct))]
    for row in rmdct: print(list(map(round, row)))

    print("qmdct  (квантованные коэффициенты)")
    qmdct = qMDCT(mdct, q_s)
    for row in qmdct: print(row)

    print("rqmdct (после обратного квантования)")
    rqmdct = rqMDCT(qmdct, q_s)
    for row in rqmdct: print(list(map(round, row)))

    print("rqrmdct (восстановленный блок)")
    # В оригинале здесь стояло range(len(rqrmdct)) — опечатка (NameError),
    # исправлено на len(rqmdct).
    rqrmdct = [[rMDCT(rqmdct, k1, k2) for k2 in range(len(rqmdct[k1]))]
               for k1 in range(len(rqmdct))]
    for row in rqrmdct: print(list(map(round, row)))

    print("absolute error (MSE-подобная метрика)")
    err = sum(
        (1 / len(rmdct)) * sum(
            (1 / len(rqr)) * (t - o) ** 2
            for t, o in zip(rqr, r)
        ) for rqr, r in zip(rqrmdct, rmdct)
    )
    print(err)

# ── тестовые данные 1D 
q_s = [1,2,4,8,16]
o_s = [48,54,80,57,8]
#d1ex(o_s,q_s)

# ── тестовые данные 2D
o_s = [
    [12,10,14,9],
    [15,12,6,20],
    [12,21,8,14],
    [15,17,17,11]
]
q_s = [
    [1,1,1,1],
    [1,1,1,1],
    [1,1,1,1],
    [1,1,1,1],
]    
q_s = [
    [1,2,9,9],
    [2,3,9,9],
    [9,9,9,9],
    [9,9,9,9]
]
#d2ex(o_s,q_s)

# ═══════════════════════════════════════════════════════════════════════════
# Цветовые преобразования RGB ↔ YCbCr (BT.601)
# ═══════════════════════════════════════════════════════════════════════════

def RGBtoYCbCr(P):
    """
    RGB → YCbCr (BT.601, как в JPEG).
    Возвращает [Y, Cb, Cr]; Cb/Cr сдвинуты на +128 → диапазон [0, 255].
    """
    R, G, B = P
    return list(map(round, [
         0.2990 * R + 0.5870 * G + 0.11140 * B,          # Y
        -0.1687 * R - 0.3313 * G + 0.5000 * B + 128,     # Cb
         0.5000 * R - 0.4187 * G - 0.0813 * B + 128,     # Cr
    ]))

def YCbCrtoRGB(P):
    """
    YCbCr → RGB. Клиппинг в [0, 255] (после округления может выйти за
    диапазон). ВНИМАНИЕ: ожидает порядок [Y, Cb, Cr].
    """
    Y, Cb, Cr = P
    res = list(map(round, [
        Y + 1.402 * (Cr - 128),
        Y - 0.34414 * (Cb - 128) - 0.71414 * (Cr - 128),
        Y + 1.772 * (Cb - 128),
    ]))
    for i in range(len(res)):
        if res[i] >= 256: res[i] = 255
        if res[i] <= -1:  res[i] = 0
    return res

# ═══════════════════════════════════════════════════════════════════════════
# Утилиты для работы с матрицами
# ═══════════════════════════════════════════════════════════════════════════
 
def ftr(p_s, func):
    """
    «For-each-to-result»: применяет func к каждому пикселю матрицы p_s.
    """
    return [[func(p_s[y][x]) for x in range(len(p_s[y]))]
            for y in range(len(p_s))]

def conv_ex(p_rgb):
    """Демо: RGB → YCbCr → RGB (показать обратимость)."""
    p_ycbcr = RGBtoYCbCr(p_rgb)
    p_rycbcr = YCbCrtoRGB(p_ycbcr)
    print(p_rgb)
    print(p_ycbcr)
    print(p_rycbcr)

p_rgb = [60,60,60]
# conv_ex(p_rgb)

# ═══════════════════════════════════════════════════════════════════════════
# Субдискретизация хромы 4:2:0
# ═══════════════════════════════════════════════════════════════════════════



#for given matrix in YCbCr pixel format output tree matrix with Y matrix Cb matrix,Cr matrix
def subsampl_4_2_0(YCbCr_2d_img):
    """
    Из H×W×3 YCbCr-изображения получает три канала:
      * Y  — H×W        (без потерь, глаз к яркости чувствителен);
      * Cb — H/2 × W/2  (усреднение по 2×2);
      * Cr — H/2 × W/2  (усреднение по 2×2).

    ВНИМАНИЕ: имена каналов вводят в заблуждение. Ожидается вход в порядке
    OpenCV YCrCb, т.е. index 1 = Cr, index 2 = Cb. Поэтому переменная
    "Cb_m" фактически хранит Cr, а "Cr_m" — Cb. Свап компенсируется в
    demo_image() при сборке пикселя. Знак '-0.1' у Cb (но не у Cr) — из
    оригинала, оставлен без изменений (влияет на округление .5).
    """
    M = YCbCr_2d_img
    H, W = len(M), len(M[0])

    Y_m = [[M[y][x][0] for x in range(W)] for y in range(H)]

    Cb_m = [[round(sum([
        M[y    ][x    ][1], M[y - 1][x    ][1],
        M[y    ][x - 1][1], M[y - 1][x - 1][1],
    ]) / 4 - 0.1) for x in range(1, W, 2)] for y in range(1, H, 2)]

    Cr_m = [[round(sum([
        M[y    ][x    ][2], M[y - 1][x    ][2],
        M[y    ][x - 1][2], M[y - 1][x - 1][2],
    ]) / 4) for x in range(1, W, 2)] for y in range(1, H, 2)]

    return [Y_m, Cb_m, Cr_m]

# ── большой тестовый YCbCr-массив (из оригинала;для
#    ручного тестирования subsampl_4_2_0) ──
o_s = [
    [
        [
            60,60,60
        ],
        [
            60,60,60
        ],
        [
            60,60,60
        ],
        [
            1,1,1
        ]
    ],
    [
        [
            1,1,1
        ],
        [
            60,1,60
        ],
        [
            40,1,40
        ],
        [
            60,1,60
        ]
    ],
    [
        [
            1,1,1
        ],
        [
            1,1,1
        ],
        [
            1,1,1
        ],
        [
            1,1,1
        ]
    ],
    [
        [
            1,1,1
        ],
        [
            1,1,1
        ],
        [
            1,1,1
        ],
        [
            1,1,1
        ]
    ],
]

# ═══════════════════════════════════════════════════════════════════════════
# Точка входа
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # demo_image()                              # требует 0.jpg и 0.png
    tsv = [48, 54, 80, 57, 8]
    tqv = [1, 2, 4, 8, 16]
    # d1ex(tsv, tqv)
    tsm = [                                      # демо-блок 4×4
        [12, 10, 14,  9],
        [15, 12,  6, 20],
        [12, 21,  8, 14],
        [15, 17, 17, 11],
    ]
    tqm = [
        [1, 2, 9, 9],
        [2, 3, 9, 9],
        [9, 9, 9, 9],
        [9, 9, 9, 9],
    ]
    #d2ex(tsm, tqm)
    pass