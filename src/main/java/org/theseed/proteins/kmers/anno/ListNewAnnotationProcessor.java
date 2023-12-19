/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import org.slf4j.LoggerFactory;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.genome.iterator.GenomeSource.Type;
import org.theseed.reports.AnnotationReporter;
import org.theseed.utils.BaseReportProcessor;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;

/**
 * This command looks at two genome sources containing identical genomes that differ only in annotation.  It will list
 * the feature ID, old annotation, old subsystem, new annotation, and new subsystem for each feature in each genome.
 * This will be a huge report on the standard output.
 *
 * The positional parameters are the name of the old-annotation genome source and the new-annotation genome source.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -o	name of output file (if not STDOUT)
 *
 * --oldType	type of genome source for the old annotations (default DIR)
 * --newType	type of genome source for the new annotations (default DIR)
 * --format		output report type (default FULL)
 *
 * @author Bruce Parrello
 *
 */
public class ListNewAnnotationProcessor extends BaseReportProcessor implements AnnotationReporter.IParms {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(ListNewAnnotationProcessor.class);
    /** old-annotation genome source */
    private GenomeSource oldGenomes;
    /** new-annotation genome source */
    private GenomeSource newGenomes;
    /** output report controller */
    private AnnotationReporter reporter;

    // COMMAND-LINE OPTIONS

    /** old-annotation genome source type */
    @Option(name = "--oldType", usage = "genome source type for old-annotation genomes")
    private GenomeSource.Type oldType;

    /** new-annotation genome source type */
    @Option(name = "--newType", usage = "genome source type for new-annotation genomes")
    private GenomeSource.Type newType;

    /** format of output report */
    @Option(name = "--format", usage = "output report format")
    private AnnotationReporter.Type outType;

    /** input file or directory for old-annotation genomes */
    @Argument(index = 0, metaVar = "oldDir", usage = "genome source for old-annotation genomes", required = true)
    private File oldDir;

    /** input file or directory for new-annotation genomes */
    @Argument(index = 1, metaVar = "newDir", usage = "genome source for new-annotation genomes", required = true)
    private File newDir;


    @Override
    protected void setReporterDefaults() {
        this.oldType = GenomeSource.Type.DIR;
        this.newType = GenomeSource.Type.DIR;
        this.outType = AnnotationReporter.Type.FULL;
    }

    @Override
    protected void validateReporterParms() throws IOException, ParseFailureException {
        // Insure the old-annotation source exists.
        if (! this.oldDir.exists())
            throw new FileNotFoundException("Old-annotation source " + this.oldDir + " is not found.");
        // Insure the new-annotation source exists.
        if (! this.newDir.exists())
            throw new FileNotFoundException("New-annotation source " + this.newDir + " is not found.");
        // Connect to the genome sources.
        this.oldGenomes = this.connect(this.oldDir, this.oldType);
        this.newGenomes = this.connect(this.newDir, this.newType);
        if (this.oldGenomes.size() != this.newGenomes.size())
            log.warn("WARNING: Genome sources are different sizes!");
        // Set up the reporter.
        this.reporter = this.outType.create();
    }

    /**
     * Connect to a genome source.
     *
     * @param dir	source directory or file name
     * @param type	source type
     *
     * @return the genome source
     *
     * @throws ParseFailureException
     * @throws IOException
     */
    private GenomeSource connect(File dir, Type type) throws IOException, ParseFailureException {
        log.info("Connecting to {} genome source {}.", type, dir);
        GenomeSource retVal = type.create(dir);
        return retVal;
    }

    @Override
    protected void runReporter(PrintWriter writer) throws Exception {
        // Write out the header.
        this.reporter.startReport(this, writer);
        // Loop through the genomes.
        int processed = 0;
        int gErrors = 0;
        int fErrors = 0;
        int fCount = 0;
        final int total = this.oldGenomes.size();
        for (Genome genome : this.oldGenomes) {
            String genomeId = genome.getId();
            processed++;
            log.info("Processing genome {} of {}: {}.", processed, total, genome);
            Genome newGenome = this.newGenomes.getGenome(genomeId);
            if (newGenome == null) {
                log.error("ERROR: Genome {} not found in new-annotation library.", genomeId);
                gErrors++;
            } else {
                for (Feature feat : genome.getFeatures()) {
                    String fid = feat.getId();
                    Feature newFeat = newGenome.getFeature(fid);
                    fCount++;
                    if (newFeat == null) {
                        log.error("ERROR: Feature {} not found in new version of {}.", fid, newGenome);
                        fErrors++;
                    } else
                        this.reporter.processFeature(feat, newFeat);
                }
            }
        }
        log.info("{} features processed in {} genomes.  {} feature errors and {} genome errors.",
                fCount, processed, fErrors, gErrors);
        this.reporter.finishReport();
    }

}
