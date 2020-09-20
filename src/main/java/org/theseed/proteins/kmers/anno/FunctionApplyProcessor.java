/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.subsystems.SubsystemProjector;
import org.theseed.utils.BaseProcessor;

/**
 * This command reads in a function-mapping file produced by the "core.utils proteins" command and applies the good function mappings
 * to the genomes in a GTO directory.  The modified genomes are output to a second directory.  The function-mapping requires deleting
 * subsystem information, after which they can optionally be projected from a projector file.
 *
 * The positional parameters are the name of the function-mapping file, the name of the input genome directory, and the name of the
 * output genome directory.  The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more progress messages on the log
 *
 * --clear		erase the output directory before proceeding
 * --project	if specified, a subsystem projector file used to project subsystems after the function update
 *
 * @author Bruce Parrello
 */
public class FunctionApplyProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FunctionApplyProcessor.class);
    /** function definition map */
    private FunctionMap funMap;
    /** function conversion map (ID to description) */
    private Map<String, String> conversionMap;
    /** subsystem projector */
    private SubsystemProjector projector;

    // COMMAND-LINE OPTIONS

    /** subsystem projector file (optional) */
    @Option(name = "--project", metaVar = "projector.tbl", usage = "if specified, a file used to project new subsystems before output")
    private File projectorFile;

    @Option(name = "--clear", usage = "clear output directory before processing")
    private boolean clearFlag;

    /** function-mapping file */
    @Argument(index = 0, metaVar = "functionMapping.tbl", usage = "function-mapping file from core.utils", required = true)
    private File conversionFile;

    /** input genome directory */
    @Argument(index = 1, metaVar = "inDir", usage = "input GTO directory", required = true)
    private File inDir;

    /** output genome directory */
    @Argument(index = 2, metaVar = "outDir", usage = "output directory", required = true)
    private File outDir;

    @Override
    protected void setDefaults() {
        this.clearFlag = false;
        this.projectorFile = null;
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Validate the input directory.
        if (! this.inDir.isDirectory())
            throw new FileNotFoundException("Input directory " + this.inDir + " is not found or invalid.");
        // Read in the mapping file.  This also creates the function map.
        this.readConversionFile(this.conversionFile);
        // Validate the projector file.
        if (this.projectorFile != null) {
            log.info("Loading subsystem projector from {}.", this.projectorFile);
            this.projector = SubsystemProjector.Load(this.projectorFile);
        }
        // Now set up the output directory.
        if (this.outDir.exists()) {
            if (! this.outDir.isDirectory())
                throw new FileNotFoundException("Output directory " + this.outDir + " is invalid.");
            if (this.clearFlag) {
                log.info("Erasing output directory {}.", this.outDir);
                FileUtils.cleanDirectory(this.outDir);
            } else
                log.info("Output will be to directory {}.", this.outDir);
        } else {
            // Here we must create the output directory.
            log.info("Creating output directory {}.", this.outDir);
            FileUtils.forceMkdir(this.outDir);
        }
        return true;
    }

    /**
     * Create the function definition map and the conversion map.
     *
     * @param conversionFileIn	file containing the conversion mapping
     *
     * @throws IOException
     */
    private void readConversionFile(File conversionFileIn) throws IOException {
        // Create the maps.
        this.funMap = new FunctionMap();
        this.conversionMap = new HashMap<String, String>(50000);
        // Read through the input file.
        try (TabbedLineReader conversionStream = new TabbedLineReader(conversionFileIn)) {
            int oldFunIdx = conversionStream.findField("patric_function");
            int newFunIdx = conversionStream.findField("core_function");
            int goodFlagIdx = conversionStream.findField("good");
            for (TabbedLineReader.Line line : conversionStream) {
                if (line.getFlag(goodFlagIdx)) {
                    // Here we have a good mapping.
                    String oldFunDesc = line.get(oldFunIdx);
                    Function oldFun = this.funMap.findOrInsert(oldFunDesc);
                    String newFunDesc = line.get(newFunIdx);
                    // Verify that this is a real function change.
                    Function newFun = this.funMap.getByName(newFunDesc);
                    if (newFun == null || ! newFun.equals(oldFun))
                        this.conversionMap.put(oldFun.getId(), newFunDesc);
                }
            }
            log.info("{} function mappings found.", this.conversionMap.size());
        }
    }

    @Override
    protected void runCommand() throws Exception {
        // Loop through the input directory, producing output genomes.
        int genomesIn = 0;
        int fInTotal = 0;
        int fChangedTotal = 0;
        log.info("Scanning input directory {}.", this.inDir);
        GenomeDirectory genomes = new GenomeDirectory(this.inDir);
        for (Genome genome : genomes) {
            genomesIn++;
            log.info("Processing genome #{}: {}.", genomesIn, genome);
            int fIn = 0;
            int fChanged = 0;
            for (Feature feat : genome.getFeatures()) {
                fIn++;
                String functionDesc = feat.getFunction();
                if (feat.getFunction() != null && ! feat.getFunction().isEmpty()) {
                    Function oldFun = this.funMap.getByName(functionDesc);
                    // Functions without mappings will generally not be in the definition map.
                    if (oldFun != null) {
                        String newFunction = this.conversionMap.get(oldFun.getId());
                        if (newFunction != null) {
                            feat.setFunction(newFunction);
                            fChanged++;
                        }
                    }
                }
            }
            log.info("{} features found and {} changed.", fIn, fChanged);
            fInTotal += fIn;
            fChangedTotal += fChanged;
            if (this.projector != null) {
                log.info("Updating subsystems in {}.", genome);
                this.projector.project(genome);
            } else {
                log.info("Deleting subsystems in {}.", genome);
                genome.clearSubsystems();
            }
            // Write the genome to the output.
            File outFile = new File(this.outDir, genome.getId() + ".gto");
            log.info("Saving genome to {}.", outFile);
            genome.update(outFile);
        }
        // All done.  Write the statistics.
        log.info("All done.  {} genomes processed, {} features analyzed, {} updated.", genomesIn, fInTotal, fChangedTotal);
    }

}
