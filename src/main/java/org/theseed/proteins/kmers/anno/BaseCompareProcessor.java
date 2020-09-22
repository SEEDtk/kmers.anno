/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.security.NoSuchAlgorithmException;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.genome.Genome;
import org.theseed.genome.compare.MatchGenomes;
import org.theseed.sequence.MD5Hex;
import org.theseed.utils.BaseProcessor;

/**
 * This is the base class for genome ORF-comparison processors.  The subclass determines the type of comparison engine
 * as well as controlling all the output.
 *
 * @author Bruce Parrello
 *
 */
public abstract class BaseCompareProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(FunctionCompareProcessor.class);
    /** map of MD5s to old genomes */
    private Map<String, File> md5GenomeMap;
    /** md5 computer */
    private MD5Hex md5Computer;

    // COMMAND-LINE OPTIONS

    /** old-genome directory */
    @Argument(index = 0, metaVar = "refDir", usage = "reference-genome directory", required = true)
    private File oldDir;

    @Override
    protected final void setDefaults() {
        this.setSubDefaults();
    }

    /**
     * Set the defaults for the subclass parameters.
     */
    protected abstract void setSubDefaults();

    @Override
    protected final boolean validateParms() throws IOException {
        if (! this.oldDir.isDirectory())
            throw new FileNotFoundException("Reference genome directory " + this.oldDir + " is not found or invalid.");
        try {
            this.validateSubParms();
            // Create the MD5 engine.
            this.md5Computer = new MD5Hex();
            // Create the map of old genomes.
            log.info("Scanning old-genome directory {}.", this.oldDir);
            this.md5GenomeMap = this.getCompareEngine().getMd5GenomeMap(this.oldDir);
            log.info("{} genomes found in {}.", this.md5GenomeMap.size(), this.oldDir);
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException(e);
        }
        return true;
    }

    /**
     * @return the comparison engine for this processor
     */
    protected abstract MatchGenomes getCompareEngine();

    /**
     * Validate the parameters for the subclass.
     */
    protected abstract void validateSubParms() throws IOException, NoSuchAlgorithmException;

    /**
     * @return the old-directory genome corresponding to this one
     *
     * @param genome	new-directory genome to check
     *
     * @throws UnsupportedEncodingException
     */
    protected File findOldGenome(Genome genome) throws UnsupportedEncodingException {
        String genomeMD5 = this.md5Computer.sequenceMD5(genome);
        File oldGenomeFile = this.md5GenomeMap.get(genomeMD5);
        return oldGenomeFile;
    }


}
