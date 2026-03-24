
------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------
              MACROECONOMICS III - DYNAMIC PROGRAMMING AND REAL BUSINESS CYCLES   
 
            Course taught by Jesús Bueren, codes developed by José Pedro Garcia 
                                           EUI 2023-2024                       
------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------



This set of codes accompanies the notes I wrote on dynamic programming. More specifically, they are meant to implement the three models which are fully developed there: 1) a deterministic neoclassical growth model (same as in PS2), 2) a stochastic neoclassical growth model (that of PS3 without the labour supply decision), and 3) a labour demand model with adjustment costs (not official material for the course).

The codes are not meant to be the quickest and most efficient ones, but rather to show you how you can use value function iteration to obtain numerical solutions for all the problems we pose in the models. You will also find copious notes all throughout the scripts which should guide you through what I am doing.

If you spot a mistake or you think the codes can be significantly improved let me know: josepedro.garcia@eui.eu

// José Garcia



-------------------------
FILES:
-------------------------

neoclassical_det.m: script, solves the deterministic neoclassical growth model (ch. 3 of the notes) with value function iteration.

neoclassical_stoch.m: script, solves the stochastic neoclassical growth model (ch. 4.1 of the notes) with value function iteration.

labour_demand.m: script, solves the labour demand model (ch. 4.2 of the notes) with value function iteration.
	
tauchen.m: function, implements the Tauchen method to discretise productivity/profitability in the models of ch. 4.
	




