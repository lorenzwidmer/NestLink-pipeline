#!/usr/bin/env nextflow

/* Define parameters */
params.outdir = "results"

/* Print pipeline info */
log.info """
        =================================
        N E S T L I N K   P I P E L I N E
        =================================
        Output dir: ${params.outdir}
        """
        .stripIndent()

/* Input channels */


/* Processes */
process sayHello {
    output:
        stdout
    """
    echo 'Hello World!'
    """
}

/* Workflow */
workflow {
    sayHello()
}
