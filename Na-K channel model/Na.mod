TITLE Na channel for HH model

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (mM) = (milli/liter)
}

NEURON {
    SUFFIX Na
    USEION na READ ena WRITE ina
    RANGE gna_bar,gna
    GLOBAL minf, hinf, mtau, htau
}

PARAMETER {
    gna_bar = 120 (mho/cm2)
}

ASSIGNED {
    v (mV)
    ina (mA/cm2)
    ena (mV)
    minf
    hinf
    mtau (ms)
    htau (ms)
    malpha mbeta halpha hbeta (/ms)
}

STATE {
    m
    h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina = (0.001)*gna_bar*m*m*m*h*(v - ena)
}


INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {
    rates(v)
    m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}

UNITSOFF

PROCEDURE rates(v (mV)) {
    TABLE minf, hinf, mtau, htau FROM -100 TO 100 WITH 200
    malpha=(0.1*(v+25))/(exp((v+25)/10)-1)
    mbeta=4*exp(v/18)

    halpha=0.07*exp(v/20)
    hbeta=1/(exp((v+30)/10)+1)
    
    mtau = 1/(malpha+mbeta)
    minf = malpha/(malpha+mbeta)
    
    htau = 1/(halpha+hbeta)
    hinf = halpha/(halpha+hbeta)
}

UNITSON