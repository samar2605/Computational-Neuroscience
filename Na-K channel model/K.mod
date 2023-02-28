TITLE K channel for HH model

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
}

NEURON {
    SUFFIX K
    USEION k READ ek WRITE ik
    RANGE gk_bar,gk
    GLOBAL ninf, ntau
}

PARAMETER {
    gk_bar = 0.036 (mho/cm2)
}

ASSIGNED {
    v (mV)
    ik (mA/cm2)
    ek (mV)
    ninf
    ntau (ms)
    nalpha nbeta (/ms)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik = gk_bar*n*n*n*n*(v - ek)
}


INITIAL {
    rates(v)
    n = ninf
}

DERIVATIVE states {
    rates(v)
    n' = (ninf - n)/ntau
}

UNITSOFF

PROCEDURE rates(v (mV)) {
    TABLE ninf, ntau FROM -100 TO 100 WITH 200
    nalpha=(0.01*(v+10))/(exp((v+10)/10)-1)
    nbeta=0.125*exp(v/80)

    ntau = 1/(nalpha+nbeta)
    ninf = nalpha/(nalpha+nbeta)
}

UNITSON