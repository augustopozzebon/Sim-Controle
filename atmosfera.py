 import numpy as np
 
def ISA(H,DT):
 # Funcao para calculo da atmosfera padrao ISA. Inclui variacao da
 # temperatura em relacao a condicao padrao.
 # Valida para altitudes entre 0 e 50.000 m. A referencia de altitude eh o
 # nivel medio do mar para a latitude de 45º. A gravidade de referencia eh
 # 9,80665 m/s^2.
 # Entradas:
 # H (m): altitude geopotencial
 # DT (ºC ou ºK): Variacao de temperatura em relacao a condicao padrao
 # Saidas (calculadas na altitude fornecida para a variacao de temperatura estipulada):
 # T (K): Temperatura
 # p (N/m^2): Pressao
 # rho (kg/m^3): Densidade
 # a (m/s): Velocidade do som
 # mu (kg/m.s): Viscosidade dinamica
 # Hp (m): Altitude pressao 
 ## Constantes padrao
     g0 = 9.80665
     R = 287.053
     gama = 1.4
     mu0 = 1.7894e-05
 ## Tabelas da atmosfera padrao
     HB = np.array([0,11000,20000,32000,47000])
     HT = np.array([11000,20000,32000,47000,50000])
     A = np.array([- 0.0065,0,0.001,0.0028,0])  
     TB = np.array([288.15,216.65,216.65,228.65,270.65])
     PB = np.array([101325,22632,5474.87,868.014,110.906])
 ## Correcao das tabelas de pressao e temperatura para a temperatura nao padrao
     tb = TB + DT
     pb=np.empty((5,1))
     pb[0] = PB[0]
     pb[1] = PB[0] * (TB[1] / TB[0]) ** (- g0 / (A[0] * R))
     pb[2] = pb[0] * np.exp(- g0 * (HB[2] - HB[1]) / (R * TB[1]))
     pb[3] = pb[1] * (TB[3] / TB[2]) ** (- g0 / (A[2] * R))
     pb[4] = pb[2] * (TB[4] / TB[3]) ** (- g0 / (A[3] * R))
 ## Calculos do modelo de atmosfera padrao
     if H < 0:
         print('Cuidado, altitude negativa. O modelo ISA nao eh valido. Resultados para H=0.')
         T = TB[0]
         p = PB[0]
         rho = p / (R * T)
         a = np.sqrt(gama * R * T)
         mu = mu0
         Hp = 0
     else:
         if H > 50000:
                 print('Cuidado, H>50.000 m. O modelo ISA nao eh valido. Resultados para H=50.000 m.')
                 T = TB[4]
                 p = pb[3] * np.exp(- g0 * (HT[4] - HB[4]) / (R * TB[4]))
                 rho = p / (R * T)
                 a = np.sqrt(gama * R * T)
                 mu = mu0 * (T / TB[0]) * ((TB[0] + 110) / (T + 110))
                 Hp = HB[4] - (R * TB[4] / g0) * np.log(p / pb[3])
         else:
             # Verifica em que camada a altitude fornecida se encontra
             i = 0
             while H - HT[i] > 0:
                 i = i + 1
             T = tb[i] + A[i] * (H - HB[i])
             if (i == 1) or (i == 4):
                 p = pb[i] * np.exp(- g0 * (H - HB[i]) / (R * tb[i]))
             else:
                 p = pb[i] * (T / tb[i]) ** (- g0 / (A[i] * R))
                 rho = p / (R * T)
                 a = np.sqrt(gama * R * T)
                 mu = mu0 * (T / TB[0]) * ((TB[0] + 110) / (T + 110))
             # Verifica em que camada da atmosfera padrao a pressao calculada se encontra
             i = 0
             while p - PB[i + 1] < 0:
                 i = i + 1
                 if i == 4:
                     break
             if (i == 1) or (i == 4):
                 Hp = HB[i] - (R * TB[i] / g0) * np.log(p / PB[i])
             else:
                 Hp = HB[i] + (TB[i] / A[i]) * ((p / PB[i]) ** (- A[i] * R / g0) - 1)
                 return T,p,rho,a,mu,Hp


T,p,rho,a,mu,Hp=ISA(5000,5)
print("Temperatura: %5.2f °C" % T)
print("Pressao: %5.2f Pa" % p)
print("Densidade: %5.2f kg/m^3" % rho)
print("Altitude pressao : %5.2f °C" % Hp)


