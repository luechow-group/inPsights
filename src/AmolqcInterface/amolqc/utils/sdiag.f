C
C
      SUBROUTINE SDIAG (A,MZA,MSA,D,MD,V,MZV,MSV,NB,IERR)
C PROGRAMM ZUR LOESUNG EINES STANDARDEIGENWERTPROBLEMS/ BESTIMMUNG DER
C EIGENWERTE D(NB) UND EIGENVEKTOREN V(NB,NB) DER MATRIX A(NB,NB) MIT
C HILFE DES FELDES E(160) FUER BIS ZU 20 BASENFUNKTIONEN. NB = ZAHL DER
C BASENFUNKTIONEN UND ZAHL DER ZU BESTIMMENDEN EIGENWERTE UND EIGENVEK-
C TOREN, DIE ANORDNUNG ERFOLFT NACH STEIGENDEN EIGENWERTEN.
C VERFAHREN NACH KLEINDIENST/FILBERG, 18. 12. 1986
      IMPLICIT NONE
      REAL*8 A, D, V
      INTEGER MZA, MSA, MD, MZV, MSV, NB, IERR
      REAL*8 B, C, E, F, G, H, P, R, S, EPS, TOL
      INTEGER N, N1, NI, I, II, J, J1, K, L, NUM
      DIMENSION A(MZA,MSA), D(MD), V(MZV,MSV), E(300)
C
      N = NB
      IERR = 0
C
      IF (MZA.LT.1.OR.MSA.LT.1.OR.MD.LT.1.OR.MZV.LT.1.OR.MSV.LT.1.OR.MZA
     &.LT.NB.OR.MSA.LT.NB.OR.MZV.LT.NB.OR.MD.LT.NB.OR.MSV.LT.NB.OR.NB.LE
     &.0.OR.NB.GT.320) THEN
         IERR = 1
         RETURN
      END IF
C
      NUM = 4
      EPS = 1.0D-30
      TOL = 1.0D-45
      IF (N .EQ. 1) GOTO 400
      DO 10 I = 1,N
      DO 10 J = 1,I
10    V(I,J) = A(I,J)
C
      DO 150 NI = 2,N
      II = N + 2 - NI
C
      DO 150 I = II,II
      L = I - 2
      H = 0.0
      G = V(I,I-1)
      IF (L) 140, 140, 20
20    DO 30 K = 1,L
30    H = H + V(I,K) ** 2
      S = H + G * G
      IF (S .GE. TOL) GOTO 50
      H = 0.D0
      GOTO 140
50    IF (H) 140, 140, 60
60    L = L + 1
      F = G
      G = SQRT(S)
      IF (F) 75, 75, 70
70    G = -G
75    H = S - F * G
      V(I,I-1) = F - G
      F = 0.0D0
C
      DO 110 J = 1,L
      V(J,I) = V(I,J) / H
      S = 0.0D0
      DO 80 K = 1,J
80    S = S + V(J,K) * V(I,K)
      J1 = J + 1
      IF (J1 .GT. L) GOTO 100
      DO 90 K = J1,L
90    S = S + V(K,J) * V(I,K)
100   E(J) = S / H
      F = F + S * V(J,I)
110   CONTINUE
      F = F / (H + H)
      DO 120 J = 1,L
120   E(J) = E(J) - F * V(I,J)
C
      DO 130 J = 1,L
      F = V(I,J)
      S = E(J)
      DO 130 K = 1,J
130   V(J,K) = V(J,K) - F * E(K) - V(I,K) * S
C
140   D(I) = H
150   E(I-1) = G
C
      D(1) = V(1,1)
      V(1,1) = 1.0D0
      DO 220 I = 2,N
      L = I - 1
      IF (D(I)) 200, 200, 170
170   DO 190 J = 1,L
      S = 0.D0
      DO 180 K = 1,L
180   S = S + V(I,K) * V(K,J)
      DO 190 K = 1,L
190   V(K,J) = V(K,J) - S * V(K,I)
200   D(I) = V(I,I)
      V(I,I) = 1.D0
210   DO 220 J = 1,L
      V(I,J) = 0.D0
      V(J,I) = 0.D0
220   CONTINUE
C
      B = 0.D0
      F = 0.D0
      E(N) = 0.D0
C
      DO 340 L = 1,N
      H = EPS * (ABS(D(L)) + ABS(E(L)))
      IF (H .GT. B) B = H
C
      DO 240 J = L,N
      IF (ABS(E(J)) .LE. B) GOTO 250
240   CONTINUE
C
250   IF (J .EQ. L) GOTO 340
C
260   P = (D(L+1) - D(L)) * 0.5D0 / E(L)
      R = SQRT(P * P + 1.0D0)
      IF (P) 270, 280, 280
270   P = P - R
      GOTO 290
280   P = P + R
290   H = D(L) - E(L) / P
      DO 300 I = L,N
300   D(I) = D(I) - H
      F = F + H
C
      P = D(J)
      C = 1.0D0
      S = 0.0D0
      J1 = J - 1
      DO 330 N1 = L,J1
      II = L + J1 - N1
      DO 330 I = II,II
      G = C * E(I)
      H = C * P
C
      IF (ABS(P) .LT. ABS(E(I))) GOTO 310
      C = E(I) / P
      R = SQRT(C * C + 1.0D0)
      E(I+1) = S * P * R
      S = C / R
      C = 1.0D0 / R
      GOTO 320
310   C = P / E(I)
      R = SQRT(C * C + 1.0D0)
      E(I+1) = S * E(I) * R
      S = 1.0D0 / R
      C = C / R
320   P = C * D(I) - S * G
      D(I+1) = H + S * (C * G + S * D(I))
      DO 330 K = 1,N
      H = V(K,I+1)
      V(K,I+1) = V(K,I) * S + H * C
      V(K,I) = V(K,I) * C - H * S
330   CONTINUE
C
      E(L) = S * P
      D(L) = C * P
      IF (ABS(E(L)) .GT. B) GOTO 260
C
340   D(L) = D(L) + F
C
      N1 = N - 1
      DO 380 I = 1,N1
      K = I
      P = D(I)
      J1 = I + 1
      DO 360 J = J1,N
      IF (D(J) .GE. P) GOTO 360
      K = J
      P = D(J)
360   CONTINUE
      IF (K .EQ. I) GOTO 380
      D(K) = D(I)
      D(I) = P
      DO 370 J = 1,N
      P = V(J,I)
      V(J,I) = V(J,K)
370   V(J,K) = P
380   CONTINUE
      RETURN
C
400   D(1) = A(1,1)
      V(1,1) = 1.D0
      RETURN
      END
