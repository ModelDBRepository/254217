
: Novel dendritic action potentials shape the computational 
: 	properties of human layer 2/3 cortical neurons( Gidon et al., 2019)
: 
: written by 
: Athanasia Papoutsi and Albert Gidon

NEURON {
	SUFFIX dCaAP
	RANGE w, vth
	RANGE K, tauA,tauB, A, B, D
	RANGE t_dCaAP
	RANGE sigma_diff,refract_period
	RANGE vrest : to normalize the K
	RANGE dcaap_count
	NONSPECIFIC_CURRENT i
}

UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
    (uS) = (microsiemens)
}

PARAMETER {
    w = 0 : by default the mechanism does not work.
    vth    (mV)   : spike threshold 
	refract_period
	tauA
	tauB
	D
	sigma_diff
	K
	vrest
	t_dCaAP 
	dcaap_count
}

ASSIGNED {
    v (mV)
   	i (nA)
}

STATE{
	A
	B
}
INITIAL {
	A=0
	B=0
: the first spike is a complication, because we only describe
: the inactivation mechanism phenomemologically, we just skip 
: the first spike	
	K=0 
:just make sure we are not in the refractory period	
	t_dCaAP = -refract_period 
	dcaap_count = 0
}

BREAKPOINT {	
	SOLVE dCaAP METHOD cnexp
	i =  -(A - B) * w * K
}

DERIVATIVE dCaAP{
	: derivative of this 
	: 		sigmoid(y) = 1 / (1 + exp(-y))
	: given y = (t-b)/tau
	: is 
	: 		(1/tau)*sigmoid(y) * (1 - sigmoid(y) )
	A' =  A * (1 - A) / tauA
	B' =  B * (1 - B) / tauB
}


BEFORE BREAKPOINT {
	: 1. activate dcaap after refractory period
	: 2. activate dcaap only if threshold was crossed
	: 3. ... and only in compartments where w is set to > 0
	if (t >= t_dCaAP + refract_period && v > vth && w > 0) {
		t_dCaAP = t
		A=0.001
		B=0	
		: K is set only in the begining of the dcaap.
		K = exp( -(v - vth) / (vth - vrest) / D ) 
		if(K > 1){
			K = 1
		}
		dcaap_count = dcaap_count + 1
		printf("Fired %.0f dCaAP at %f ms with K = %f\n",dcaap_count, t_dCaAP,K)
	}
	: B is equal to 0 only when spike is initially triggered.
	if(B == 0 && t >= t_dCaAP + sigma_diff && w > 0){
		B=0.001
	}		
}
