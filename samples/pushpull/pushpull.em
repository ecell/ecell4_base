@{MAIN_STEPPER='DE1'}
@{MAIN_PROCESS='MassActionFluxProcess'}

@#{MAIN_STEPPER='NR1'}
@#{MAIN_PROCESS='GillespieProcess'}

#Stepper DiscreteEventStepper( NR1 )
#{
#}

Stepper ODEStepper( DE1 )
{
}

@{
VOL = 1e-16
N_A = 6.02214e+23

def C2N( conc ):
    num = N_A * VOL * conc
    print round( num )

koff_ratio = .1

Keq = 0.03

kon = 0.02e9

#{C2N( 249e-9 )};
NS = 200
S_conc = NS / N_A / VOL

kcatkoff = Keq * S_conc * kon
kcat = kcatkoff * (1-koff_ratio)
koff = kcatkoff * koff_ratio

NSp = 0
NP = 10
NK = 10
}

System System( / )
{
        StepperID       @MAIN_STEPPER;

        Variable Variable( SIZE )
        {
                Value   @(VOL);
        }


        Variable Variable( S )
        {
                Value   @(NS);
        }

        Variable Variable( P )
        {
                Value   @(NP);
        }

        Variable Variable( K )
        {
                Value   @(NK);
        }

        Variable Variable( KS )
        {
                Value   0;
        }

        Variable Variable( Sp )
        {
                Value   @(NSp);
        }

        Variable Variable( PSp )
        {
                Value   0;
        }

        Process @(MAIN_PROCESS)( R1 )
        {
                VariableReferenceList   [ _ :.:S      -1 ] 
                                        [ _ :.:K  -1 ]
                                        [ _ :.:KS  1];
                k       @(kon);
        }

        Process @(MAIN_PROCESS)( R2 )
        {
                VariableReferenceList   [ _ :.:KS -1 ]
                                        [ _ :.:S       1 ] 
                                        [ _ :.:K   1 ];
                k       @(koff);
        }

        Process @(MAIN_PROCESS)( R3 )
        {
                VariableReferenceList   [ _ :.:KS -1 ]
                                        [ _ :.:Sp      1 ] 
                                        [ _ :.:K   1 ];
                k       @(kcat);
        }


        Process @(MAIN_PROCESS)( R4 )
        {
                VariableReferenceList   [ _ :.:Sp       -1 ]
                                        [ _ :.:P     -1 ] 
                                        [ _ :.:PSp   1 ];
                k       @(kon);
        }

        Process @(MAIN_PROCESS)( R5 )
        {
                VariableReferenceList   [ _ :.:PSp  -1 ]
                                        [ _ :.:Sp         1 ]
                                        [ _ :.:P      1 ];
                k       @(koff);
        }

        Process @(MAIN_PROCESS)( R6 )
        {
                VariableReferenceList   [ _ :.:PSp  -1 ]
                                        [ _ :.:P        1 ]
                                        [ _ :.:S      1 ]; 
                k       @(kcat);
        }


        

        

}

