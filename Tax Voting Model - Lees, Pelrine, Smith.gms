$TITLE  Final-Project: 3 Consumer, 2 Industry Voting Model

SETS I /A,B,C/;

PARAMETERS
  ALPHA (I)  share of X in consumer utility (Cobb-Douglas)
  SIZE (I)   proportional group size
  LENDOW (I) endowment of labor
  KENDOW(I)  endowment of capital
  BETA      share of K in production of X
  GAMMA       share of K in production of Y
  CP (I)     Level of consumer preference (0<CP<1) ;

ALPHA("A") = .9;
ALPHA("B") = .1;
ALPHA("C") = .2;
SIZE("A") = .24;
SIZE("B") = .48;
SIZE("C") = .28;
LENDOW("A") = 0;
LENDOW("B") = 10;
LENDOW("C") = 9;
KENDOW("A") = 10;
KENDOW("B") = 0;
KENDOW("C") = 1;
BETA = .25;
GAMMA = .75;

VARIABLES
  TOT_WELF total welfare
  WA,WB,WC welfare of respective consumer for use as objective variable;

POSITIVE VARIABLES
  TAX tax on good X
  X activity level for sector X
  Y activity level for sector Y
  K activity level for capital market
  L activity level for labor market
  PX price of X
  PY price of Y
  PK price of capital
  PL price of labor
  CONS(I) income of consumer I
  W(I) welfare of consumer I;

EQUATIONS

  PRF_X zero-profit for sector X
  PRF_Y zero-profit for sector Y

  MKT_X market clearing condition for commodity X
  MKT_Y market clearing condition for commodity Y
  MKT_K market clearing condition for factor K
  MKT_L market clearing condition for factor L

  I_CONS(I) income for consumer I

  WELFARE(I) welfare of consumer I
  TOT_WELFARE total welfare
  WELFA,WELFB,WELFC welfare objective functions;

* ZERO PROFIT CONDITIONS
  PRF_X..        (1 + TAX) * (PK*(((BETA*PL)/((1-BETA)*PK))**(1-BETA)) + PL*((((1-BETA)*PK)/(BETA*PL))**BETA))
                 =G= PX;
  PRF_Y..        PK*(((GAMMA*PL)/((1-GAMMA)*PK))**(1-GAMMA)) + PL*((((1-GAMMA)*PK)/(GAMMA*PL))**GAMMA)
                 =G= PY;

* MARKET CLEARING CONDITIONS
  MKT_X..        X =G= SUM(I, SIZE(I)*CONS(I)*ALPHA(I)/PX);
  MKT_Y..        Y =G= SUM(I, SIZE(I)*CONS(I)*(1 - ALPHA(I))/PY);
  MKT_K..        SUM(I, KENDOW(I)*SIZE(I)) =G= (X*((BETA*PL)/((1-BETA)*PK))**(1-BETA)) + (Y*((GAMMA*PL)/((1-GAMMA)*PK))**(1-GAMMA));
  MKT_L..        SUM(I, LENDOW(I)*SIZE(I)) =G= (X*(((1-BETA)*PK)/(BETA*PL))**BETA) + (Y*(((1-GAMMA)*PK)/(GAMMA*PL))**GAMMA);

  I_CONS(I)..     CONS(I) =E= (LENDOW(I)*PL + KENDOW(I)*PK + SIZE(I)*TAX*X*PX/(1+TAX));

  WELFARE(I)..    W(I) =E= (((CONS(I)*ALPHA(I)/PX)**ALPHA(I)) * ((CONS(I)*(1 - ALPHA(I))/PY)**(1-ALPHA(I))));
  TOT_WELFARE..   TOT_WELF =E= PROD(I, W(I)**SIZE(I));
  WELFA.. WA =E= W("A");
  WELFB.. WB =E= W("B");
  WELFC.. WC =E= W("C");


MODEL EQUIL /PRF_X.X, PRF_Y.Y,
            MKT_X.PX, MKT_Y.PY, MKT_K.PK, MKT_L.PL,
            I_CONS.CONS,WELFARE.W,TOT_WELFARE,WELFA,WELFB,WELFC/;

X.L = 1.5;
Y.L = 1.5;
PX.L = 2;
PY.L = 2;
PK.L = 1;
CONS.L(I) = 1;
W.L(I) = 1;
TAX.L = .1;
PL.FX = 1;

OPTION MPEC = nlpec;
SOLVE EQUIL USING MPEC maximizing TOT_WELF;

* VOTER TURNOUT PARAMETER CALCULATION

SETS H /A,B,C/;

PARAMETERS
 WW(H,I)
 WELFGAP(I)
 MINGAP
 MAXGAP
 VOTING(I);

SOLVE EQUIL USING MPEC maximizing WA;

WW("A","A") = W.L("A")/SIZE("A"); WW("A","B") = W.L("B")/SIZE("B"); WW("A","C") = W.L("C")/SIZE("C");

SOLVE EQUIL USING MPEC maximizing WB;

WW("B","A") = W.L("A")/SIZE("A"); WW("B","B") = W.L("B")/SIZE("B"); WW("B","C") = W.L("C")/SIZE("C");

SOLVE EQUIL USING MPEC maximizing WC;

WW("C","A") = W.L("A")/SIZE("A"); WW("C","B") = W.L("B")/SIZE("B"); WW("C","C") = W.L("C")/SIZE("C");

If((WW("A","A") - WW("B","A")) >= (WW("A","A") - WW("C","A")),
WELFGAP("A") = (WW("A","A") - WW("B","A"));
ELSE
WELFGAP("A") = (WW("A","A") - WW("C","A"));)

If((WW("B","B") - WW("A","B")) >= (WW("B","B") - WW("C","B")),
WELFGAP("B") = (WW("B","B") - WW("A","B"));
ELSE
WELFGAP("B") = (WW("B","B") - WW("C","B"));)

If((WW("C","C") - WW("A","C")) >= abs(WW("C","C") - WW("B","C")),
WELFGAP("C") = (WW("C","C") - WW("A","C"));
ELSE
WELFGAP("C") = (WW("C","C") - WW("B","C"));)

if(smin(I,WELFGAP(I)) = WELFGAP("A"),
MINGAP = WELFGAP("A");
ELSE IF (smin(I,WELFGAP(I)) = WELFGAP("B"),
MINGAP = WELFGAP("B");
ELSE
MINGAP = WELFGAP("C");));

if(smax(I,WELFGAP(I)) = WELFGAP("A"),
MAXGAP = WELFGAP("A");
ELSE IF (smax(I,WELFGAP(I)) = WELFGAP("B"),
MAXGAP = WELFGAP("B");
ELSE
MAXGAP = WELFGAP("C");));

VOTING(I) = 1/(1.3 + 2.71828**((-4) * (WELFGAP(I) - MINGAP) / (MAXGAP - MINGAP) + 0.9)) ;

*Determines Condorcet loser, with 0 = none, 1 = A, 2 = B, 3 = C.

SETS P /NOVOTE,VOTE/;

PARAMETER Condorcet_Loser(P);
Condorcet_Loser("NOVOTE") = 0;
Condorcet_Loser("NOVOTE")$(SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) > SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) AND
         SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")) > SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Loser("NOVOTE")$(SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) > SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) AND
         SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) > SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Loser("NOVOTE")$(SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) > SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) AND
         SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) > SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

Condorcet_Loser("VOTE") = 0;
Condorcet_Loser("VOTE")$(VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) > VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) AND
         VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")) > VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Loser("VOTE")$(VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) > VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) AND
         VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) > VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Loser("VOTE")$(VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) > VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) AND
         VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) > VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

*Determines Condorcet winner, with 0 = none, 1 = A, 2 = B, 3 = C.

PARAMETER Condorcet_Winner(P);
Condorcet_Winner("NOVOTE") = 0;
Condorcet_Winner("NOVOTE")$(SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) < SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) AND
         SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")) < SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Winner("NOVOTE")$(SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) < SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) AND
         SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) < SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Winner("NOVOTE")$(SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) < SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) AND
         SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) < SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

Condorcet_Winner("VOTE") = 0;
Condorcet_Winner("VOTE")$(VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) < VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) AND
         VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")) < VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Winner("VOTE")$(VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) < VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) AND
         VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) < VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Winner("VOTE")$(VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) < VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) AND
         VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) < VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;


SETS V /NOVOTE,FPPNOTURN,FPPTURN,RONOTURN,ROTURN/

PARAMETERS
 ACTUAL_WINNER(V)
 RESULTS1(*,V);

ACTUAL_WINNER("NOVOTE") = 0;

*First-past-the-post: Chooses largest group, solves for their optimal level.
*In relation to a real-world situation, this assumes each group puts forward their optimal tax as an option in the voting.

* NO VOTING ALLOWED - BASE CASE

SOLVE EQUIL USING MPEC maximizing TOT_WELF;

RESULTS1("ACTUAL_WINNER","NOVOTE") = ACTUAL_WINNER("NOVOTE");
RESULTS1("CONDORCET_WINNER","NOVOTE") = 0;
RESULTS1("CONDORCET_LOSER","NOVOTE") = 0;
RESULTS1("VOTENUMA","NOVOTE") = 0;
RESULTS1("VOTENUMB","NOVOTE") = 0;
RESULTS1("VOTENUMC","NOVOTE") = 0;
RESULTS1("TOT_WELF","NOVOTE") = TOT_WELF.L;
RESULTS1("WA","NOVOTE") = WA.L;
RESULTS1("WB","NOVOTE") = WB.L;
RESULTS1("WC","NOVOTE") = WC.L;
RESULTS1("TAX","NOVOTE") = TAX.L;
RESULTS1("X","NOVOTE") = X.L;
RESULTS1("Y","NOVOTE") = Y.L;
RESULTS1("PX","NOVOTE") = PX.L;
RESULTS1("PY","NOVOTE") = PY.L;
RESULTS1("PK","NOVOTE") = PK.L;
RESULTS1("PL","NOVOTE") = PL.L;

* WITHOUT VOTER TURNOUT PARAMETER

If(smax(I,SIZE(I)) = SIZE("A"),
ACTUAL_WINNER("FPPNOTURN") = 1;
SOLVE EQUIL USING MPEC maximizing WA;
ELSE If(smax(I,SIZE(I)) = SIZE("B"),
ACTUAL_WINNER("FPPNOTURN") = 2;
SOLVE EQUIL USING MPEC maximizing WB;
ELSE
ACTUAL_WINNER("FPPNOTURN") = 3;
SOLVE EQUIL USING MPEC maximizing WC;));

RESULTS1("ACTUAL_WINNER","FPPNOTURN") = ACTUAL_WINNER("FPPNOTURN");
RESULTS1("CONDORCET_WINNER","FPPNOTURN") = CONDORCET_WINNER("NOVOTE");
RESULTS1("CONDORCET_LOSER","FPPNOTURN") = CONDORCET_LOSER("NOVOTE");
RESULTS1("VOTENUMA","FPPNOTURN") = SIZE("A");
RESULTS1("VOTENUMB","FPPNOTURN") = SIZE("B");
RESULTS1("VOTENUMC","FPPNOTURN") = SIZE("C");
RESULTS1("TOT_WELF","FPPNOTURN") = TOT_WELF.L;
RESULTS1("WA","FPPNOTURN") = WA.L;
RESULTS1("WB","FPPNOTURN") = WB.L;
RESULTS1("WC","FPPNOTURN") = WC.L;
RESULTS1("TAX","FPPNOTURN") = TAX.L;
RESULTS1("X","FPPNOTURN") = X.L;
RESULTS1("Y","FPPNOTURN") = Y.L;
RESULTS1("PX","FPPNOTURN") = PX.L;
RESULTS1("PY","FPPNOTURN") = PY.L;
RESULTS1("PK","FPPNOTURN") = PK.L;
RESULTS1("PL","FPPNOTURN") = PL.L;

* WITH VOTER TURNOUT PARAMETER

If(smax(I,SIZE(I)*VOTING(I)) = SIZE("A")*VOTING("A"),
ACTUAL_WINNER("FPPTURN") = 1;
SOLVE EQUIL USING MPEC maximizing WA;
ELSE If(smax(I,SIZE(I)*VOTING(I)) = SIZE("B")*VOTING("B"),
ACTUAL_WINNER("FPPTURN") = 2;
SOLVE EQUIL USING MPEC maximizing WB;
ELSE
ACTUAL_WINNER("FPPTURN") = 3;
SOLVE EQUIL USING MPEC maximizing WC;));

RESULTS1("ACTUAL_WINNER","FPPTURN") = ACTUAL_WINNER("FPPTURN");
RESULTS1("CONDORCET_WINNER","FPPTURN") = CONDORCET_WINNER("VOTE");
RESULTS1("CONDORCET_LOSER","FPPTURN") = CONDORCET_LOSER("VOTE");
RESULTS1("VOTENUMA","FPPTURN") = VOTING("A")*SIZE("A");
RESULTS1("VOTENUMB","FPPTURN") = VOTING("B")*SIZE("B");
RESULTS1("VOTENUMC","FPPTURN") = VOTING("C")*SIZE("C");
RESULTS1("TOT_WELF","FPPTURN") = TOT_WELF.L;
RESULTS1("WA","FPPTURN") = WA.L;
RESULTS1("WB","FPPTURN") = WB.L;
RESULTS1("WC","FPPTURN") = WC.L;
RESULTS1("TAX","FPPTURN") = TAX.L;
RESULTS1("X","FPPTURN") = X.L;
RESULTS1("Y","FPPTURN") = Y.L;
RESULTS1("PX","FPPTURN") = PX.L;
RESULTS1("PY","FPPTURN") = PY.L;
RESULTS1("PK","FPPTURN") = PK.L;
RESULTS1("PL","FPPTURN") = PL.L;

*Runoff Voting: In the first round, the group with the least votes is eliminated.
*They then decide which of the two remaining options they prefer, and add their votes (in this case, their size) to that group.

* WITHOUT VOTER TURNOUT PARAMETER

If (smin(I,SIZE(I)) = SIZE("A"),
RESULTS1("VOTENUMA","RONOTURN") = 0;
RESULTS1("VOTENUMB","RONOTURN") = SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"));
RESULTS1("VOTENUMC","RONOTURN") = SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A"));
         If (SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) > SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")),
         ACTUAL_WINNER("RONOTURN") = 2;
         SOLVE EQUIL USING MPEC maximizing WB;
         ELSE
         ACTUAL_WINNER("RONOTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE IF (smin(I,SIZE(I)) = SIZE("B"),
RESULTS1("VOTENUMB","RONOTURN") = 0;
RESULTS1("VOTENUMA","RONOTURN") = SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"));
RESULTS1("VOTENUMC","RONOTURN") = SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"));
         If (SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) > SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")),
         ACTUAL_WINNER("RONOTURN") = 1;
         SOLVE EQUIL USING MPEC maximizing WA;
         ELSE
         ACTUAL_WINNER("RONOTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE
RESULTS1("VOTENUMC","RONOTURN") = 0;
RESULTS1("VOTENUMA","RONOTURN") = SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C"));
RESULTS1("VOTENUMB","RONOTURN") = SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C"));
        If (SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) > SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")),
        ACTUAL_WINNER("RONOTURN") = 1;
        SOLVE EQUIL USING MPEC maximizing WA;
        ELSE
        ACTUAL_WINNER("RONOTURN") = 2;
        SOLVE EQUIL USING MPEC maximizing WB;)
););

RESULTS1("ACTUAL_WINNER","RONOTURN") = ACTUAL_WINNER("RONOTURN");
RESULTS1("CONDORCET_WINNER","RONOTURN") = CONDORCET_WINNER("NOVOTE");
RESULTS1("CONDORCET_LOSER","RONOTURN") = CONDORCET_LOSER("NOVOTE");
RESULTS1("TOT_WELF","RONOTURN") = TOT_WELF.L;
RESULTS1("WA","RONOTURN") = WA.L;
RESULTS1("WB","RONOTURN") = WB.L;
RESULTS1("WC","RONOTURN") = WC.L;
RESULTS1("TAX","RONOTURN") = TAX.L;
RESULTS1("X","RONOTURN") = X.L;
RESULTS1("Y","RONOTURN") = Y.L;
RESULTS1("PX","RONOTURN") = PX.L;
RESULTS1("PY","RONOTURN") = PY.L;
RESULTS1("PK","RONOTURN") = PK.L;
RESULTS1("PL","RONOTURN") = PL.L;

* WITH VOTER TURNOUT PARAMETER

If (smin(I,VOTING(I)*SIZE(I)) = VOTING("A")*SIZE("A"),
RESULTS1("VOTENUMA","ROTURN") = 0;
RESULTS1("VOTENUMB","ROTURN") = VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"));
RESULTS1("VOTENUMC","ROTURN") = VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A"));
         If (VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) > VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")),
         ACTUAL_WINNER("ROTURN") = 2;
         SOLVE EQUIL USING MPEC maximizing WB;
         ELSE
         ACTUAL_WINNER("ROTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE IF (smin(I,VOTING(I)*SIZE(I)) = VOTING("B")*SIZE("B"),
RESULTS1("VOTENUMB","ROTURN") = 0;
RESULTS1("VOTENUMA","ROTURN") = VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"));
RESULTS1("VOTENUMC","ROTURN") = VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"));
         If (VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) > VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")),
         ACTUAL_WINNER("ROTURN") = 1;
         SOLVE EQUIL USING MPEC maximizing WA;
         ELSE
         ACTUAL_WINNER("ROTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE
RESULTS1("VOTENUMC","ROTURN") = 0;
RESULTS1("VOTENUMA","ROTURN") = VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C"));
RESULTS1("VOTENUMB","ROTURN") = VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C"));
        If (VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) > VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")),
        ACTUAL_WINNER("ROTURN") = 1;
        SOLVE EQUIL USING MPEC maximizing WA;
        ELSE
        ACTUAL_WINNER("ROTURN") = 2;
        SOLVE EQUIL USING MPEC maximizing WB;)
););

RESULTS1("ACTUAL_WINNER","ROTURN") = ACTUAL_WINNER("ROTURN");
RESULTS1("CONDORCET_WINNER","ROTURN") = CONDORCET_WINNER("VOTE");
RESULTS1("CONDORCET_LOSER","ROTURN") = CONDORCET_LOSER("VOTE");
RESULTS1("TOT_WELF","ROTURN") = TOT_WELF.L;
RESULTS1("WA","ROTURN") = WA.L;
RESULTS1("WB","ROTURN") = WB.L;
RESULTS1("WC","ROTURN") = WC.L;
RESULTS1("TAX","ROTURN") = TAX.L;
RESULTS1("X","ROTURN") = X.L;
RESULTS1("Y","ROTURN") = Y.L;
RESULTS1("PX","ROTURN") = PX.L;
RESULTS1("PY","ROTURN") = PY.L;
RESULTS1("PK","ROTURN") = PK.L;
RESULTS1("PL","ROTURN") = PL.L;

DISPLAY RESULTS1;

*Initial values show FPP: Condorcet winner loses and loser wins, RO: Condorcet winner loses.
*Counterfactual: change sizes so loser wins/winner loses on FPP with turnout:

SIZE("A") = .21;
SIZE("B") = .6;
SIZE("C") = .19;

SOLVE EQUIL USING MPEC maximizing TOT_WELF;

* VOTER TURNOUT PARAMETER CALCULATION


SOLVE EQUIL USING MPEC maximizing WA;

WW("A","A") = W.L("A")/SIZE("A"); WW("A","B") = W.L("B")/SIZE("B"); WW("A","C") = W.L("C")/SIZE("C");

SOLVE EQUIL USING MPEC maximizing WB;

WW("B","A") = W.L("A")/SIZE("A"); WW("B","B") = W.L("B")/SIZE("B"); WW("B","C") = W.L("C")/SIZE("C");

SOLVE EQUIL USING MPEC maximizing WC;

WW("C","A") = W.L("A")/SIZE("A"); WW("C","B") = W.L("B")/SIZE("B"); WW("C","C") = W.L("C")/SIZE("C");

If((WW("A","A") - WW("B","A")) >= (WW("A","A") - WW("C","A")),
WELFGAP("A") = (WW("A","A") - WW("B","A"));
ELSE
WELFGAP("A") = (WW("A","A") - WW("C","A"));)

If((WW("B","B") - WW("A","B")) >= (WW("B","B") - WW("C","B")),
WELFGAP("B") = (WW("B","B") - WW("A","B"));
ELSE
WELFGAP("B") = (WW("B","B") - WW("C","B"));)

If((WW("C","C") - WW("A","C")) >= abs(WW("C","C") - WW("B","C")),
WELFGAP("C") = (WW("C","C") - WW("A","C"));
ELSE
WELFGAP("C") = (WW("C","C") - WW("B","C"));)

if(smin(I,WELFGAP(I)) = WELFGAP("A"),
MINGAP = WELFGAP("A");
ELSE IF (smin(I,WELFGAP(I)) = WELFGAP("B"),
MINGAP = WELFGAP("B");
ELSE
MINGAP = WELFGAP("C");));

if(smax(I,WELFGAP(I)) = WELFGAP("A"),
MAXGAP = WELFGAP("A");
ELSE IF (smax(I,WELFGAP(I)) = WELFGAP("B"),
MAXGAP = WELFGAP("B");
ELSE
MAXGAP = WELFGAP("C");));

VOTING(I) = 1/(1.3 + 2.71828**((-4) * (WELFGAP(I) - MINGAP) / (MAXGAP - MINGAP) + 0.9)) ;

PARAMETER Condorcet_Loser(P);
Condorcet_Loser("NOVOTE") = 0;
Condorcet_Loser("NOVOTE")$(SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) > SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) AND
         SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")) > SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Loser("NOVOTE")$(SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) > SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) AND
         SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) > SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Loser("NOVOTE")$(SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) > SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) AND
         SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) > SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

Condorcet_Loser("VOTE") = 0;
Condorcet_Loser("VOTE")$(VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) > VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) AND
         VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")) > VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Loser("VOTE")$(VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) > VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) AND
         VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) > VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Loser("VOTE")$(VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) > VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) AND
         VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) > VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

PARAMETER Condorcet_Winner(P);
Condorcet_Winner("NOVOTE") = 0;
Condorcet_Winner("NOVOTE")$(SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) < SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) AND
         SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")) < SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Winner("NOVOTE")$(SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) < SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) AND
         SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) < SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Winner("NOVOTE")$(SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) < SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) AND
         SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) < SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

Condorcet_Winner("VOTE") = 0;
Condorcet_Winner("VOTE")$(VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) < VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) AND
         VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")) < VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Winner("VOTE")$(VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) < VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) AND
         VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) < VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Winner("VOTE")$(VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) < VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) AND
         VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) < VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;



PARAMETERS
 RESULTS2(*,V);

ACTUAL_WINNER("NOVOTE") = 0;

* NO VOTING ALLOWED - BASE CASE

SOLVE EQUIL USING MPEC maximizing TOT_WELF;

RESULTS2("ACTUAL_WINNER","NOVOTE") = ACTUAL_WINNER("NOVOTE");
RESULTS2("CONDORCET_WINNER","NOVOTE") = 0;
RESULTS2("CONDORCET_LOSER","NOVOTE") = 0;
RESULTS2("VOTENUMA","NOVOTE") = 0;
RESULTS2("VOTENUMB","NOVOTE") = 0;
RESULTS2("VOTENUMC","NOVOTE") = 0;
RESULTS2("TOT_WELF","NOVOTE") = TOT_WELF.L;
RESULTS2("WA","NOVOTE") = WA.L;
RESULTS2("WB","NOVOTE") = WB.L;
RESULTS2("WC","NOVOTE") = WC.L;
RESULTS2("TAX","NOVOTE") = TAX.L;
RESULTS2("X","NOVOTE") = X.L;
RESULTS2("Y","NOVOTE") = Y.L;
RESULTS2("PX","NOVOTE") = PX.L;
RESULTS2("PY","NOVOTE") = PY.L;
RESULTS2("PK","NOVOTE") = PK.L;
RESULTS2("PL","NOVOTE") = PL.L;

* WITHOUT VOTER TURNOUT PARAMETER

If(smax(I,SIZE(I)) = SIZE("A"),
ACTUAL_WINNER("FPPNOTURN") = 1;
SOLVE EQUIL USING MPEC maximizing WA;
ELSE If(smax(I,SIZE(I)) = SIZE("B"),
ACTUAL_WINNER("FPPNOTURN") = 2;
SOLVE EQUIL USING MPEC maximizing WB;
ELSE
ACTUAL_WINNER("FPPNOTURN") = 3;
SOLVE EQUIL USING MPEC maximizing WC;));

RESULTS2("ACTUAL_WINNER","FPPNOTURN") = ACTUAL_WINNER("FPPNOTURN");
RESULTS2("CONDORCET_WINNER","FPPNOTURN") = CONDORCET_WINNER("NOVOTE");
RESULTS2("CONDORCET_LOSER","FPPNOTURN") = CONDORCET_LOSER("NOVOTE");
RESULTS2("VOTENUMA","FPPNOTURN") = SIZE("A");
RESULTS2("VOTENUMB","FPPNOTURN") = SIZE("B");
RESULTS2("VOTENUMC","FPPNOTURN") = SIZE("C");
RESULTS2("TOT_WELF","FPPNOTURN") = TOT_WELF.L;
RESULTS2("WA","FPPNOTURN") = WA.L;
RESULTS2("WB","FPPNOTURN") = WB.L;
RESULTS2("WC","FPPNOTURN") = WC.L;
RESULTS2("TAX","FPPNOTURN") = TAX.L;
RESULTS2("X","FPPNOTURN") = X.L;
RESULTS2("Y","FPPNOTURN") = Y.L;
RESULTS2("PX","FPPNOTURN") = PX.L;
RESULTS2("PY","FPPNOTURN") = PY.L;
RESULTS2("PK","FPPNOTURN") = PK.L;
RESULTS2("PL","FPPNOTURN") = PL.L;

* WITH VOTER TURNOUT PARAMETER

If(smax(I,SIZE(I)*VOTING(I)) = SIZE("A")*VOTING("A"),
ACTUAL_WINNER("FPPTURN") = 1;
SOLVE EQUIL USING MPEC maximizing WA;
ELSE If(smax(I,SIZE(I)*VOTING(I)) = SIZE("B")*VOTING("B"),
ACTUAL_WINNER("FPPTURN") = 2;
SOLVE EQUIL USING MPEC maximizing WB;
ELSE
ACTUAL_WINNER("FPPTURN") = 3;
SOLVE EQUIL USING MPEC maximizing WC;));

RESULTS2("ACTUAL_WINNER","FPPTURN") = ACTUAL_WINNER("FPPTURN");
RESULTS2("CONDORCET_WINNER","FPPTURN") = CONDORCET_WINNER("VOTE");
RESULTS2("CONDORCET_LOSER","FPPTURN") = CONDORCET_LOSER("VOTE");
RESULTS2("VOTENUMA","FPPTURN") = VOTING("A")*SIZE("A");
RESULTS2("VOTENUMB","FPPTURN") = VOTING("B")*SIZE("B");
RESULTS2("VOTENUMC","FPPTURN") = VOTING("C")*SIZE("C");
RESULTS2("TOT_WELF","FPPTURN") = TOT_WELF.L;
RESULTS2("WA","FPPTURN") = WA.L;
RESULTS2("WB","FPPTURN") = WB.L;
RESULTS2("WC","FPPTURN") = WC.L;
RESULTS2("TAX","FPPTURN") = TAX.L;
RESULTS2("X","FPPTURN") = X.L;
RESULTS2("Y","FPPTURN") = Y.L;
RESULTS2("PX","FPPTURN") = PX.L;
RESULTS2("PY","FPPTURN") = PY.L;
RESULTS2("PK","FPPTURN") = PK.L;
RESULTS2("PL","FPPTURN") = PL.L;

* WITHOUT VOTER TURNOUT PARAMETER

If (smin(I,SIZE(I)) = SIZE("A"),
RESULTS2("VOTENUMA","RONOTURN") = 0;
RESULTS2("VOTENUMB","RONOTURN") = SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"));
RESULTS2("VOTENUMC","RONOTURN") = SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A"));
         If (SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) > SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")),
         ACTUAL_WINNER("RONOTURN") = 2;
         SOLVE EQUIL USING MPEC maximizing WB;
         ELSE
         ACTUAL_WINNER("RONOTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE IF (smin(I,SIZE(I)) = SIZE("B"),
RESULTS2("VOTENUMB","RONOTURN") = 0;
RESULTS2("VOTENUMA","RONOTURN") = SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"));
RESULTS2("VOTENUMC","RONOTURN") = SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"));
         If (SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) > SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")),
         ACTUAL_WINNER("RONOTURN") = 1;
         SOLVE EQUIL USING MPEC maximizing WA;
         ELSE
         ACTUAL_WINNER("RONOTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE
RESULTS2("VOTENUMC","RONOTURN") = 0;
RESULTS2("VOTENUMA","RONOTURN") = SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C"));
RESULTS2("VOTENUMB","RONOTURN") = SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C"));
        If (SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) > SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")),
        ACTUAL_WINNER("RONOTURN") = 1;
        SOLVE EQUIL USING MPEC maximizing WA;
        ELSE
        ACTUAL_WINNER("RONOTURN") = 2;
        SOLVE EQUIL USING MPEC maximizing WB;)
););

RESULTS2("ACTUAL_WINNER","RONOTURN") = ACTUAL_WINNER("RONOTURN");
RESULTS2("CONDORCET_WINNER","RONOTURN") = CONDORCET_WINNER("NOVOTE");
RESULTS2("CONDORCET_LOSER","RONOTURN") = CONDORCET_LOSER("NOVOTE");
RESULTS2("TOT_WELF","RONOTURN") = TOT_WELF.L;
RESULTS2("WA","RONOTURN") = WA.L;
RESULTS2("WB","RONOTURN") = WB.L;
RESULTS2("WC","RONOTURN") = WC.L;
RESULTS2("TAX","RONOTURN") = TAX.L;
RESULTS2("X","RONOTURN") = X.L;
RESULTS2("Y","RONOTURN") = Y.L;
RESULTS2("PX","RONOTURN") = PX.L;
RESULTS2("PY","RONOTURN") = PY.L;
RESULTS2("PK","RONOTURN") = PK.L;
RESULTS2("PL","RONOTURN") = PL.L;

* WITH VOTER TURNOUT PARAMETER

If (smin(I,VOTING(I)*SIZE(I)) = VOTING("A")*SIZE("A"),
RESULTS2("VOTENUMA","ROTURN") = 0;
RESULTS2("VOTENUMB","ROTURN") = VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"));
RESULTS2("VOTENUMC","ROTURN") = VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A"));
         If (VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) > VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")),
         ACTUAL_WINNER("ROTURN") = 2;
         SOLVE EQUIL USING MPEC maximizing WB;
         ELSE
         ACTUAL_WINNER("ROTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE IF (smin(I,VOTING(I)*SIZE(I)) = VOTING("B")*SIZE("B"),
RESULTS2("VOTENUMB","ROTURN") = 0;
RESULTS2("VOTENUMA","ROTURN") = VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"));
RESULTS2("VOTENUMC","ROTURN") = VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"));
         If (VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) > VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")),
         ACTUAL_WINNER("ROTURN") = 1;
         SOLVE EQUIL USING MPEC maximizing WA;
         ELSE
         ACTUAL_WINNER("ROTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE
RESULTS2("VOTENUMC","ROTURN") = 0;
RESULTS2("VOTENUMA","ROTURN") = VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C"));
RESULTS2("VOTENUMB","ROTURN") = VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C"));
        If (VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) > VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")),
        ACTUAL_WINNER("ROTURN") = 1;
        SOLVE EQUIL USING MPEC maximizing WA;
        ELSE
        ACTUAL_WINNER("ROTURN") = 2;
        SOLVE EQUIL USING MPEC maximizing WB;)
););

RESULTS2("ACTUAL_WINNER","ROTURN") = ACTUAL_WINNER("ROTURN");
RESULTS2("CONDORCET_WINNER","ROTURN") = CONDORCET_WINNER("VOTE");
RESULTS2("CONDORCET_LOSER","ROTURN") = CONDORCET_LOSER("VOTE");
RESULTS2("TOT_WELF","ROTURN") = TOT_WELF.L;
RESULTS2("WA","ROTURN") = WA.L;
RESULTS2("WB","ROTURN") = WB.L;
RESULTS2("WC","ROTURN") = WC.L;
RESULTS2("TAX","ROTURN") = TAX.L;
RESULTS2("X","ROTURN") = X.L;
RESULTS2("Y","ROTURN") = Y.L;
RESULTS2("PX","ROTURN") = PX.L;
RESULTS2("PY","ROTURN") = PY.L;
RESULTS2("PK","ROTURN") = PK.L;
RESULTS2("PL","ROTURN") = PL.L;


DISPLAY RESULTS2;

*Counterfactual: change sizes so Condorcet winner loses on RO with turnout:

SIZE("A") = .20;
SIZE("B") = .56;
SIZE("C") = .24;

SOLVE EQUIL USING MPEC maximizing TOT_WELF;

* VOTER TURNOUT PARAMETER CALCULATION

SOLVE EQUIL USING MPEC maximizing WA;

WW("A","A") = W.L("A")/SIZE("A"); WW("A","B") = W.L("B")/SIZE("B"); WW("A","C") = W.L("C")/SIZE("C");

SOLVE EQUIL USING MPEC maximizing WB;

WW("B","A") = W.L("A")/SIZE("A"); WW("B","B") = W.L("B")/SIZE("B"); WW("B","C") = W.L("C")/SIZE("C");

SOLVE EQUIL USING MPEC maximizing WC;

WW("C","A") = W.L("A")/SIZE("A"); WW("C","B") = W.L("B")/SIZE("B"); WW("C","C") = W.L("C")/SIZE("C");

If((WW("A","A") - WW("B","A")) >= (WW("A","A") - WW("C","A")),
WELFGAP("A") = (WW("A","A") - WW("B","A"));
ELSE
WELFGAP("A") = (WW("A","A") - WW("C","A"));)

If((WW("B","B") - WW("A","B")) >= (WW("B","B") - WW("C","B")),
WELFGAP("B") = (WW("B","B") - WW("A","B"));
ELSE
WELFGAP("B") = (WW("B","B") - WW("C","B"));)

If((WW("C","C") - WW("A","C")) >= abs(WW("C","C") - WW("B","C")),
WELFGAP("C") = (WW("C","C") - WW("A","C"));
ELSE
WELFGAP("C") = (WW("C","C") - WW("B","C"));)

if(smin(I,WELFGAP(I)) = WELFGAP("A"),
MINGAP = WELFGAP("A");
ELSE IF (smin(I,WELFGAP(I)) = WELFGAP("B"),
MINGAP = WELFGAP("B");
ELSE
MINGAP = WELFGAP("C");));

if(smax(I,WELFGAP(I)) = WELFGAP("A"),
MAXGAP = WELFGAP("A");
ELSE IF (smax(I,WELFGAP(I)) = WELFGAP("B"),
MAXGAP = WELFGAP("B");
ELSE
MAXGAP = WELFGAP("C");));

VOTING(I) = 1/(1.3 + 2.71828**((-4) * (WELFGAP(I) - MINGAP) / (MAXGAP - MINGAP) + 0.9)) ;

PARAMETER Condorcet_Loser(P);
Condorcet_Loser("NOVOTE") = 0;
Condorcet_Loser("NOVOTE")$(SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) > SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) AND
         SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")) > SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Loser("NOVOTE")$(SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) > SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) AND
         SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) > SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Loser("NOVOTE")$(SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) > SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) AND
         SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) > SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

Condorcet_Loser("VOTE") = 0;
Condorcet_Loser("VOTE")$(VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) > VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) AND
         VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")) > VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Loser("VOTE")$(VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) > VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) AND
         VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) > VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Loser("VOTE")$(VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) > VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) AND
         VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) > VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

PARAMETER Condorcet_Winner(P);
Condorcet_Winner("NOVOTE") = 0;
Condorcet_Winner("NOVOTE")$(SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) < SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) AND
         SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")) < SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Winner("NOVOTE")$(SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) < SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")) AND
         SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) < SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Winner("NOVOTE")$(SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) < SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")) AND
         SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) < SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;

Condorcet_Winner("VOTE") = 0;
Condorcet_Winner("VOTE")$(VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) < VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) AND
         VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")) < VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"))) = 1;
Condorcet_Winner("VOTE")$(VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) < VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")) AND
         VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) < VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"))) = 2;
Condorcet_Winner("VOTE")$(VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) < VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")) AND
         VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) < VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"))) = 3;



PARAMETERS
 RESULTS3(*,V);

ACTUAL_WINNER("NOVOTE") = 0;

* NO VOTING ALLOWED - BASE CASE

SOLVE EQUIL USING MPEC maximizing TOT_WELF;

RESULTS3("ACTUAL_WINNER","NOVOTE") = ACTUAL_WINNER("NOVOTE");
RESULTS3("CONDORCET_WINNER","NOVOTE") = 0;
RESULTS3("CONDORCET_LOSER","NOVOTE") = 0;
RESULTS3("VOTENUMA","NOVOTE") = 0;
RESULTS3("VOTENUMB","NOVOTE") = 0;
RESULTS3("VOTENUMC","NOVOTE") = 0;
RESULTS3("TOT_WELF","NOVOTE") = TOT_WELF.L;
RESULTS3("WA","NOVOTE") = WA.L;
RESULTS3("WB","NOVOTE") = WB.L;
RESULTS3("WC","NOVOTE") = WC.L;
RESULTS3("TAX","NOVOTE") = TAX.L;
RESULTS3("X","NOVOTE") = X.L;
RESULTS3("Y","NOVOTE") = Y.L;
RESULTS3("PX","NOVOTE") = PX.L;
RESULTS3("PY","NOVOTE") = PY.L;
RESULTS3("PK","NOVOTE") = PK.L;
RESULTS3("PL","NOVOTE") = PL.L;

* WITHOUT VOTER TURNOUT PARAMETER

If(smax(I,SIZE(I)) = SIZE("A"),
ACTUAL_WINNER("FPPNOTURN") = 1;
SOLVE EQUIL USING MPEC maximizing WA;
ELSE If(smax(I,SIZE(I)) = SIZE("B"),
ACTUAL_WINNER("FPPNOTURN") = 2;
SOLVE EQUIL USING MPEC maximizing WB;
ELSE
ACTUAL_WINNER("FPPNOTURN") = 3;
SOLVE EQUIL USING MPEC maximizing WC;));

RESULTS3("ACTUAL_WINNER","FPPNOTURN") = ACTUAL_WINNER("FPPNOTURN");
RESULTS3("CONDORCET_WINNER","FPPNOTURN") = CONDORCET_WINNER("NOVOTE");
RESULTS3("CONDORCET_LOSER","FPPNOTURN") = CONDORCET_LOSER("NOVOTE");
RESULTS3("VOTENUMA","FPPNOTURN") = SIZE("A");
RESULTS3("VOTENUMB","FPPNOTURN") = SIZE("B");
RESULTS3("VOTENUMC","FPPNOTURN") = SIZE("C");
RESULTS3("TOT_WELF","FPPNOTURN") = TOT_WELF.L;
RESULTS3("WA","FPPNOTURN") = WA.L;
RESULTS3("WB","FPPNOTURN") = WB.L;
RESULTS3("WC","FPPNOTURN") = WC.L;
RESULTS3("TAX","FPPNOTURN") = TAX.L;
RESULTS3("X","FPPNOTURN") = X.L;
RESULTS3("Y","FPPNOTURN") = Y.L;
RESULTS3("PX","FPPNOTURN") = PX.L;
RESULTS3("PY","FPPNOTURN") = PY.L;
RESULTS3("PK","FPPNOTURN") = PK.L;
RESULTS3("PL","FPPNOTURN") = PL.L;

* WITH VOTER TURNOUT PARAMETER

If(smax(I,SIZE(I)*VOTING(I)) = SIZE("A")*VOTING("A"),
ACTUAL_WINNER("FPPTURN") = 1;
SOLVE EQUIL USING MPEC maximizing WA;
ELSE If(smax(I,SIZE(I)*VOTING(I)) = SIZE("B")*VOTING("B"),
ACTUAL_WINNER("FPPTURN") = 2;
SOLVE EQUIL USING MPEC maximizing WB;
ELSE
ACTUAL_WINNER("FPPTURN") = 3;
SOLVE EQUIL USING MPEC maximizing WC;));

RESULTS3("ACTUAL_WINNER","FPPTURN") = ACTUAL_WINNER("FPPTURN");
RESULTS3("CONDORCET_WINNER","FPPTURN") = CONDORCET_WINNER("VOTE");
RESULTS3("CONDORCET_LOSER","FPPTURN") = CONDORCET_LOSER("VOTE");
RESULTS3("VOTENUMA","FPPTURN") = VOTING("A")*SIZE("A");
RESULTS3("VOTENUMB","FPPTURN") = VOTING("B")*SIZE("B");
RESULTS3("VOTENUMC","FPPTURN") = VOTING("C")*SIZE("C");
RESULTS3("TOT_WELF","FPPTURN") = TOT_WELF.L;
RESULTS3("WA","FPPTURN") = WA.L;
RESULTS3("WB","FPPTURN") = WB.L;
RESULTS3("WC","FPPTURN") = WC.L;
RESULTS3("TAX","FPPTURN") = TAX.L;
RESULTS3("X","FPPTURN") = X.L;
RESULTS3("Y","FPPTURN") = Y.L;
RESULTS3("PX","FPPTURN") = PX.L;
RESULTS3("PY","FPPTURN") = PY.L;
RESULTS3("PK","FPPTURN") = PK.L;
RESULTS3("PL","FPPTURN") = PL.L;

* WITHOUT VOTER TURNOUT PARAMETER

If (smin(I,SIZE(I)) = SIZE("A"),
RESULTS3("VOTENUMA","RONOTURN") = 0;
RESULTS3("VOTENUMB","RONOTURN") = SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A"));
RESULTS3("VOTENUMC","RONOTURN") = SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A"));
         If (SIZE("B") + SIZE("A")$(WW("B","A") > WW("C","A")) > SIZE("C") + SIZE("A")$(WW("C","A") > WW("B","A")),
         ACTUAL_WINNER("RONOTURN") = 2;
         SOLVE EQUIL USING MPEC maximizing WB;
         ELSE
         ACTUAL_WINNER("RONOTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE IF (smin(I,SIZE(I)) = SIZE("B"),
RESULTS3("VOTENUMB","RONOTURN") = 0;
RESULTS3("VOTENUMA","RONOTURN") = SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B"));
RESULTS3("VOTENUMC","RONOTURN") = SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B"));
         If (SIZE("A") + SIZE("B")$(WW("A","B") > WW("C","B")) > SIZE("C") + SIZE("B")$(WW("C","B") > WW("A","B")),
         ACTUAL_WINNER("RONOTURN") = 1;
         SOLVE EQUIL USING MPEC maximizing WA;
         ELSE
         ACTUAL_WINNER("RONOTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE
RESULTS3("VOTENUMC","RONOTURN") = 0;
RESULTS3("VOTENUMA","RONOTURN") = SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C"));
RESULTS3("VOTENUMB","RONOTURN") = SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C"));
        If (SIZE("A") + SIZE("C")$(WW("A","C") > WW("B","C")) > SIZE("B") + SIZE("C")$(WW("B","C") > WW("A","C")),
        ACTUAL_WINNER("RONOTURN") = 1;
        SOLVE EQUIL USING MPEC maximizing WA;
        ELSE
        ACTUAL_WINNER("RONOTURN") = 2;
        SOLVE EQUIL USING MPEC maximizing WB;)
););

RESULTS3("ACTUAL_WINNER","RONOTURN") = ACTUAL_WINNER("RONOTURN");
RESULTS3("CONDORCET_WINNER","RONOTURN") = CONDORCET_WINNER("NOVOTE");
RESULTS3("CONDORCET_LOSER","RONOTURN") = CONDORCET_LOSER("NOVOTE");
RESULTS3("TOT_WELF","RONOTURN") = TOT_WELF.L;
RESULTS3("WA","RONOTURN") = WA.L;
RESULTS3("WB","RONOTURN") = WB.L;
RESULTS3("WC","RONOTURN") = WC.L;
RESULTS3("TAX","RONOTURN") = TAX.L;
RESULTS3("X","RONOTURN") = X.L;
RESULTS3("Y","RONOTURN") = Y.L;
RESULTS3("PX","RONOTURN") = PX.L;
RESULTS3("PY","RONOTURN") = PY.L;
RESULTS3("PK","RONOTURN") = PK.L;
RESULTS3("PL","RONOTURN") = PL.L;

* WITH VOTER TURNOUT PARAMETER

If (smin(I,VOTING(I)*SIZE(I)) = VOTING("A")*SIZE("A"),
RESULTS3("VOTENUMA","ROTURN") = 0;
RESULTS3("VOTENUMB","ROTURN") = VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A"));
RESULTS3("VOTENUMC","ROTURN") = VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A"));
         If (VOTING("B")*SIZE("B") + VOTING("A")*SIZE("A")$(WW("B","A") > WW("C","A")) > VOTING("C")*SIZE("C") + VOTING("A")*SIZE("A")$(WW("C","A") > WW("B","A")),
         ACTUAL_WINNER("ROTURN") = 2;
         SOLVE EQUIL USING MPEC maximizing WB;
         ELSE
         ACTUAL_WINNER("ROTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE IF (smin(I,VOTING(I)*SIZE(I)) = VOTING("B")*SIZE("B"),
RESULTS3("VOTENUMB","ROTURN") = 0;
RESULTS3("VOTENUMA","ROTURN") = VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B"));
RESULTS3("VOTENUMC","ROTURN") = VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B"));
         If (VOTING("A")*SIZE("A") + VOTING("B")*SIZE("B")$(WW("A","B") > WW("C","B")) > VOTING("C")*SIZE("C") + VOTING("B")*SIZE("B")$(WW("C","B") > WW("A","B")),
         ACTUAL_WINNER("ROTURN") = 1;
         SOLVE EQUIL USING MPEC maximizing WA;
         ELSE
         ACTUAL_WINNER("ROTURN") = 3;
         SOLVE EQUIL USING MPEC maximizing WC;)
ELSE
RESULTS3("VOTENUMC","ROTURN") = 0;
RESULTS3("VOTENUMA","ROTURN") = VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C"));
RESULTS3("VOTENUMB","ROTURN") = VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C"));
        If (VOTING("A")*SIZE("A") + VOTING("C")*SIZE("C")$(WW("A","C") > WW("B","C")) > VOTING("B")*SIZE("B") + VOTING("C")*SIZE("C")$(WW("B","C") > WW("A","C")),
        ACTUAL_WINNER("ROTURN") = 1;
        SOLVE EQUIL USING MPEC maximizing WA;
        ELSE
        ACTUAL_WINNER("ROTURN") = 2;
        SOLVE EQUIL USING MPEC maximizing WB;)
););

RESULTS3("ACTUAL_WINNER","ROTURN") = ACTUAL_WINNER("ROTURN");
RESULTS3("CONDORCET_WINNER","ROTURN") = CONDORCET_WINNER("VOTE");
RESULTS3("CONDORCET_LOSER","ROTURN") = CONDORCET_LOSER("VOTE");
RESULTS3("TOT_WELF","ROTURN") = TOT_WELF.L;
RESULTS3("WA","ROTURN") = WA.L;
RESULTS3("WB","ROTURN") = WB.L;
RESULTS3("WC","ROTURN") = WC.L;
RESULTS3("TAX","ROTURN") = TAX.L;
RESULTS3("X","ROTURN") = X.L;
RESULTS3("Y","ROTURN") = Y.L;
RESULTS3("PX","ROTURN") = PX.L;
RESULTS3("PY","ROTURN") = PY.L;
RESULTS3("PK","ROTURN") = PK.L;
RESULTS3("PL","ROTURN") = PL.L;

DISPLAY RESULTS1;
DISPLAY RESULTS2;
DISPLAY RESULTS3;
