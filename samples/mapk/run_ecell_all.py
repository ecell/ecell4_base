import os

#MODEL_FILE = 'model.eml'
ESS_FILE = 'run_ecell_model4.py'

# Register jobs.

aJobIDList = []

for i in xrange(0,5):
        
        VALUE_OF_S = i * 1000
        aParameterDict = { 'MODEL_FILE': MODEL_FILE, 'VALUE_OF_S': VALUE_OF_S }

        #registerEcellSession( ESS file, parameters, files that ESS uses )
        aJobID = registerEcellSession( ESS_FILE, aParameterDict, [ MODEL_FILE, ])
        aJobIDList.append( aJobID ) # Memorize the job IDs in aJobIDList.

# Run the registered jobs.

run()

for aJobID in aJobIDList: 

        print " --- job id = %s ---" % aJobID
        print getStdout( aJobID )  # Print the output of each job. 
