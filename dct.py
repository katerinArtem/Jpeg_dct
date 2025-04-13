
from math import cos,pi,degrees,sqrt,acos

def cust(N,i,k):return cos((pi/N)*(i+0.5)*k)
def DCT(k,g_s):return sum([g_s[i]*cust(len(g_s),i,k) for i in range(len(g_s))])
def rDCT(u,G_s):return (2/len(G_s))*(G_s[0]*(1/2) + sum([G_s[i]*cos((pi/len(G_s)*(u+0.5)*i)) for i in range(1,len(G_s))]))
def QDCT(k,G_s,q_s):return round(G_s[k]/q_s[k])
def rQDCT(k,t_s,q_s):return t_s[k]*q_s[k]

def d1ex(o_s,q_s):
    print("S:",o_s)
    G_s = [DCT(m,o_s) for m in range(0,len(o_s))]
    print("D:",G_s)
    print('Q:',q_s)
    t_s = [QDCT(k,G_s,q_s) for k in range(0,len(o_s))]
    print('T:',t_s)
    rt_s = [rQDCT(k,t_s,q_s) for k in range(0,len(o_s))]
    print('D\':',rt_s)
    rtG_s = [rDCT(u,rt_s) for u in range(0,len(o_s))]
    print("S\'",list(map(round,rtG_s)))
    rG_s = [rDCT(u,G_s) for u in range(0,len(o_s))]
    print("S",list(map(round,rG_s)))

def MDCT(o_s,k1,k2):return sum([DCT(k2,o_s[i1])*cust(len(o_s),i1,k1) for i1 in range(len(o_s))])
def rMDCT(G_s,k1,k2):return  (2/len(G_s))*(rDCT(k2,G_s[0])*(1/2) + sum([rDCT(k2,G_s[i1])*cos((pi/len(G_s))*(k1+0.5)*i1) for i1 in range(1,len(G_s))]))
def qMDCT(G_s,Q_s):return [[round(g/q) for g,q in zip(G,Q)]for G,Q in zip(G_s,Q_s)]
def rqMDCT(T_s,Q_s):return [[g*q for g,q in zip(G,Q)]for G,Q in zip(T_s,Q_s)]

def d2ex(o_s,q_s):
    print('orig')
    for i in o_s:print(i)

    print("mdct")
    mdct = [[MDCT(o_s,k1,k2) for k2 in range(len(o_s[k1]))] for k1 in range(len(o_s))]
    for i in mdct:print(i)

    print("rmdct")
    rmdct = [[rMDCT(mdct,k1,k2) for k2 in range(len(mdct[k1]))]  for k1 in range(len(mdct))]
    for i in rmdct:print(list(map(round,i)))
            
    print("qmdct")
    qmdct = qMDCT(mdct,q_s) 
    for i in qmdct:print(i)


    print("rqmdct")
    rqmdct = rqMDCT(qmdct,q_s) 
    for i in qmdct:print(list(map(round,i)))


    print("rqrmdct")
    rqrmdct = [[rMDCT(rqmdct,k1,k2) for k2 in range(len(rqmdct[k1]))]  for k1 in range(len(rqmdct))] 
    for i in rqrmdct:print(list(map(round,i)))

    print("absolute error")
    print(sum([(1/len(rmdct))*sum([ ((1/len(rqr))*(t-o)**2) for t,o in zip(rqr,r)]) for rqr,r in zip(rqrmdct,rmdct)]))

q_s = [1,2,4,8,16]
o_s = [48,54,80,57,8]

#d1ex(o_s,q_s)

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


def RGBtoYCbCr(P:list):
    R,G,B = P
    return list(map(round,[
            0.2990*R + 0.5870*G + 0.11140*B,  
            -0.1687*R - 0.3313*G + 0.5000*B + 128, 
            0.5000*R - 0.4187*G - 0.0813*B + 128,    
    ]))

def YCbCrtoRGB(P:list):
    Y,Cb,Cr = P
    return list(map(round,[
             Y + 1.402*(Cr-128),  
             Y - 0.34414*(Cb-128) - 0.71414*(Cr-128), 
             Y + 1.772*(Cb-128),  
    ]))
 
def ftr(p_s,func):return [[func(p_s[y][x]) for x in range(len(p_s[y]))]for y in range(len(p_s))] 

def conv_ex(p_rgb):
    p_ycbcr = RGBtoYCbCr(p_rgb)
    p_rycbcr = YCbCrtoRGB(p_ycbcr)
    print(p_rgb)
    print(p_ycbcr)
    print(p_rycbcr)

p_rgb = [60,60,60]

#for given matrix in YCbCr pixel format output tree matrix with Y matrix Cb matrix,Cr matrix
def subsampl_4_2_0(YCbCr_2d_img):
    M = YCbCr_2d_img
    Cb_m,Cr_m = [],[]
    Y_m = [[s[0] for s in l]for l in M]
    Cb_m = [[round(sum([
        M[y][x][1],
        M[y-1][x][1],
        M[y][x-1][1],
        M[y-1][x-1][1]
        ])/4) for x in range(1,len(M[y]),2)] for y in range(1,len(M),2)]
    Cr_m = [[round(sum([
        M[y][x][2],
        M[y-1][x][2],
        M[y][x-1][2],
        M[y-1][x-1][2]
        ])/4) for x in range(1,len(M[y]),2)] for y in range(1,len(M),2)]
    return [Y_m,Cb_m,Cr_m]

#conv_ex(p_rgb)

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

