TITLE  K-A channel


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

NEURON {
	SUFFIX ka
	USEION k READ ek WRITE ik
        RANGE gkabar,gka
        GLOBAL minf,hinf,tauh,taum
}

PARAMETER {
	gkabar=.01 (mho/cm2)
        vhalfm=-33.6   (mV)
        vhalfh=-83   (mV)
        a0l=0.08      (/ms)
        a0n=0.02    (/ms)
        zetam=-3    (1)
        zetah=4    (1)
        gmm=0.6   (1)
        gmh=1   (1)
}

STATE {
	m
        h
}

INITIAL {
        rates(v)
        m=minf
        h=hinf
}

ASSIGNED {
        v (mV)
	ik (mA/cm2)
        ek (mV)
	celsius (degC)
        minf
        hinf      
        tauh (ms)
        taum (ms)
        gka (mho/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gkabar*m*h
	ik = gka*(v-ek)

}


FUNCTION alpm(v(mV)) {
  alpm = exp(1.e-3*zetam*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betm(v(mV)) {
  betm = exp(1.e-3*zetam*gmm*(v-vhalfm)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION alph(v(mV)) {
  alph = exp(1.e-3*zetah*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION beth(v(mV)) {
  beth = exp(1.e-3*zetah*gmh*(v-vhalfh)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states { 
        rates(v)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,q10
        q10=3^((celsius-30)/10)
        a = alpm(v)
        minf = 1/(1 + a)
        taum = betm(v)/(q10*a0n*(1+a))
        a = alph(v)
        hinf = 1/(1+ a)
        tauh = beth(v)/(q10*a0l*(1 + a))
}







