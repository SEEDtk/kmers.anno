package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileFilter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.basic.BaseProcessor;
import org.theseed.basic.ParseFailureException;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.SubsystemRow;
import org.theseed.genome.iterator.GenomeSource;
import org.theseed.io.FieldInputStream;
import org.theseed.io.FieldInputStream.Record;
import org.theseed.io.JsonListInputStream;

import com.github.cliftonlabs.json_simple.JsonArray;
import com.github.cliftonlabs.json_simple.JsonObject;

/**
 * This sub-command will update the BV-BRC JSON dump files with new annotations in GTOs for the
 * same genomes.  This is a highly specialized command.  For each genome sub-directory with the
 * JSON dump data in it, some files will be copied verbatim, some will require minor field
 * substitutions, and some will need to be completely rebuilt.  Our basic strategy will be to
 * copy the verbatim files first.  Then we will loop through genome_feature.json, replacing the
 * "product" field for each feature.  In parallel, we will rebuild an entry in the subsystem.json
 * file for each feature with subsystems.
 *
 * The positional parameters are the name of the master directory for the JSON genome dumps,
 * the name of the input genome source, and the name of the master directory for the output
 * JSON genome dumps.
 *
 * The command-line options are as follows.
 *
 * -h	display command-line usage
 * -v	display more frequent log messages
 * -t	type of input genome source (default DIR)
 *
 * --clear	erase the output directory before processing
 */
public class UpdateJsonProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(UpdateJsonProcessor.class);
    /** input genome source */
    private GenomeSource genomes;
    /** input JSON directories */
    private File[] genomeDirs;
    /** list of genome_feature field names, in column order */
    private String[] featureCols;
    /** map of genome_feature field names to types */
    private static final Map<String, JsonType> FEATURE_FIELDS = Map.ofEntries(
                    Map.entry("patric_id", JsonType.STRING),
                    Map.entry("public",  JsonType.BOOLEAN),
                    Map.entry("genome_name", JsonType.STRING),
                    Map.entry("genome_id", JsonType.STRING),
                    Map.entry("product", JsonType.STRING),
                    Map.entry("feature_type", JsonType.STRING),
                    Map.entry("accession", JsonType.STRING),
                    Map.entry("strand", JsonType.STRING),
                    Map.entry("start", JsonType.INTEGER),
                    Map.entry("end", JsonType.INTEGER),
                    Map.entry("location", JsonType.STRING),
                    Map.entry("aa_sequence_md5", JsonType.STRING),
                    Map.entry("aa_length", JsonType.INTEGER),
                    Map.entry("na_sequence_md5", JsonType.STRING),
                    Map.entry("na_length", JsonType.INTEGER),
                    Map.entry("refseq_locus_tag", JsonType.STRING),
                    Map.entry("gene", JsonType.STRING),
                    Map.entry("gene_id", JsonType.STRING),
                    Map.entry("annotation", JsonType.STRING),
                    Map.entry("protein_id", JsonType.STRING),
                    Map.entry("segments", JsonType.LIST),
                    Map.entry("taxon_id", JsonType.INTEGER)
                );
    /** set of file names for straight copy */
    private static final Set<String> COPY_FILES = Set.of("genome.json", "protein_structure.json",
            "sp_gene.json");
    /** genome ID string pattern */
    private static final Pattern GENOME_PATTERN = Pattern.compile("\\d+\\.\\d+");
    /** file filter for subdirectories of the input JSON dump that are genome IDs */
    protected static final FileFilter JGENOME_DIR_FILTER = new FileFilter() {

        @Override
        public boolean accept(File pathname) {
            // Insure we are a directory.
            boolean retVal = pathname.isDirectory();
            if (retVal) {
                String dirName = pathname.getName();
                retVal = GENOME_PATTERN.matcher(dirName).matches();
            }
            return retVal;
        }

    };

    // COMMAND-LINE OPTIONS

    /** if specified, the output directory will be erased before processingg */
    @Option(name = "--clear", usage = "erase the output directory before processing")
    private boolean clearFlag;

    /** input genome source type */
    @Option(name = "--type", aliases = {"-t" }, usage = "input genome source type")
    private GenomeSource.Type sourceType;

    /** input JSON master directory */
    @Argument(index = 0, metaVar = "jsonInDir", usage = "JSON dump input master directory", required = true)
    private File jsonInDir;

    /** input genome source */
    @Argument(index = 1, metaVar = "genomeInDir", usage = "input genome source with new annotation data", required = true)
    private File genomeInDir;

    /** output directort */
    @Argument(index = 2, metaVar = "jsonOutDir", usage = "JSON dump output master directory", required = true)
    private File jsonOutDir;

    /**
     * This enum describes a JSON field type.
     */
    protected static enum JsonType {
        BOOLEAN {
            @Override
            public Object valueOf(Record record, int idx) {
                return record.getFlag(idx);
            }
        }, INTEGER {
            @Override
            public Object valueOf(Record record, int idx) {
                return record.getInt(idx);
            }
        }, FLOAT {
            @Override
            public Object valueOf(Record record, int idx) {
                return record.getDouble(idx);
            }
        }, STRING {
            @Override
            public Object valueOf(Record record, int idx) {
                return record.get(idx);
            }
        }, LIST {
            @Override
            public Object valueOf(Record record, int idx) {
                return record.getList(idx);
            }
        };

        /**
         * Convert a JSON record field of this type into an output field.
         *
         * @param obj		output JSON object
         * @param name		field name
         * @param record	record containing the field
         * @param idx		index of the field in the record
         */
        public void output(JsonObject obj, String name, FieldInputStream.Record record, int idx) {
            obj.put(name, this.valueOf(record, idx));
        }

        /**
         * Format a field for storage in a JSON object.
         *
         * @param record	record containing the field
         * @param idx		index of the field iin the record
         *
         * @return the object to put in the JSON field
         */
        public abstract Object valueOf(FieldInputStream.Record record, int idx);

    }

    @Override
    protected void setDefaults() {
        this.clearFlag = false;
        this.sourceType = GenomeSource.Type.DIR;
    }

    @Override
    protected boolean validateParms() throws IOException, ParseFailureException {
        // Insure the input directory is valid.
        if (! this.jsonInDir.isDirectory())
            throw new FileNotFoundException("Input JSON directory " + this.jsonInDir
                    + " is not found or invalid.");
        log.info("Scanning input directory {}.", this.jsonInDir);
        this.genomeDirs = this.jsonInDir.listFiles(JGENOME_DIR_FILTER);
        if (this.genomeDirs.length <= 0)
            throw new IOException("No genome subdirectories found in " + this.jsonInDir + ".");
        // Get the genome source and verify the genomes in it.
        log.info("Connecting to {} genome source in {},", this.sourceType, this.genomeInDir);
        this.genomes = this.sourceType.create(this.genomeInDir);
        // Insure all the input genomes are in the source.
        int badGenomes = 0;
        Set<String> genomeIds = this.genomes.getIDs();
        for (File gDir : this.genomeDirs) {
            String genomeId = gDir.getName();
            if (! genomeIds.contains(genomeId)) {
                log.error("Missing genome {} in {}.", genomeId, this.genomeInDir);
                badGenomes++;
            }
        }
        if (badGenomes > 0)
            throw new ParseFailureException(Integer.toString(badGenomes) + " genomes from "
                    + this.jsonInDir + " not found in " + this.genomeInDir + ".");
        // Set up the output directory.
        if (! this.jsonOutDir.isDirectory()) {
            log.info("Creating output directory {}.", this.jsonOutDir);
            FileUtils.forceMkdir(this.jsonOutDir);
        } else if (this.clearFlag) {
            log.info("Erasing output directory {}.", this.jsonOutDir);
            FileUtils.cleanDirectory(this.jsonOutDir);
        } else
            log.info("Output wll be to {},", this.jsonOutDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        final int nGenomes = this.genomeDirs.length;
        log.info("{} genomes to process.", nGenomes);
        // Set up some counters.
        int gCount = 0;
        int substitutions = 0;
        int subRecords = 0;
        int copies = 0;
        // Loop through the input genome JSON dump directories.
        for (File genomeDir : this.genomeDirs) {
            // Get the GTO.
            String genomeId = genomeDir.getName();
            Genome genome = this.genomes.getGenome(genomeId);
            gCount++;
            log.info("Processing genome {} of {}: {}.", gCount, nGenomes, genome);
            // Create the output directory.
            File outDir = new File(this.jsonOutDir, genomeId);
            if (! outDir.isDirectory()) {
                log.info("Creating output directory {}.", outDir);
                FileUtils.forceMkdir(outDir);
            }
            // Copy all the verbatim files.
            for (String name : COPY_FILES) {
                File inFile = new File(genomeDir, name);
                File outFile = new File(outDir, name);
                if (inFile.exists()) {
                    FileUtils.copyFile(inFile, outFile);
                    copies++;
                }
            }
            // Get the input feature file name, and the output feature and subsystem
            // file names.
            File inFile =  new File(genomeDir, "genome_feature.json");
            File outFile = new File(outDir, "genome_feature.json");
            File subFile = new File(outDir, "subsystem.json");
            // Set up the output subsystem JSON object.
            JsonArray subArray = new JsonArray();
            try (JsonListInputStream inStream = new JsonListInputStream(inFile)) {
                // Set up the input fields.
                this.featureCols = new String[FEATURE_FIELDS.size()];
                for (var featureField : FEATURE_FIELDS.keySet()) {
                    int idx = inStream.findField(featureField);
                    this.featureCols[idx] = featureField;
                }
                // Two of them are special.
                int idIdx = inStream.findField("patric_id");
                int prodIdx = inStream.findField("product");
                long lastMsg = System.currentTimeMillis();
                // Start the feature output stream.
                JsonArray featArray = new JsonArray();
                while (inStream.hasNext()) {
                    FieldInputStream.Record record = inStream.next();
                    String fid = record.get(idIdx);
                    if (! StringUtils.isBlank(fid)) {
                        String product = record.get(prodIdx);
                        // Find the feature in the GTO.
                        Feature feat = genome.getFeature(fid);
                        if (feat == null)
                            log.warn("{} not found in {}.", fid, genome);
                        else {
                            // We have a feature in the GTO.  Check the function.
                            String function = feat.getPegFunction();
                            if (! function.equals(product)) {
                                // Here the function has changed.
                                record.setField("product", function);
                                substitutions++;
                            }
                            // Now check for subsystems.  For each one, we output a subsystem record.
                            var subs = feat.getSubsystemRows();
                            for (SubsystemRow sub : subs) {
                                JsonObject subObject = new JsonObject();
                                subObject.put("patric_id", fid);
                                subObject.put("role_name", this.computeRole(sub, function));
                                subObject.put("active", sub.isActive() ? "active" : "inactive");
                                subObject.put("subsystem_name", sub.getName());
                                List<String> classes = sub.getClassifications();
                                if (classes.size() >= 1)
                                    subObject.put("superclass", classes.get(0));
                                if (classes.size() >= 2)
                                    subObject.put("class", classes.get(1));
                                if (classes.size() >= 3)
                                    subObject.put("subclass", classes.get(2));
                                subArray.add(subObject);
                                subRecords++;
                            }
                        }
                    }
                    // Output the feature record.
                    JsonObject featObject = new JsonObject();
                    for (int i = 0; i < this.featureCols.length; i++) {
                        String name = this.featureCols[i];
                        JsonType type = FEATURE_FIELDS.get(name);
                        type.output(featObject, name, record, i);
                    }
                    featArray.add(featObject);
                    long current = System.currentTimeMillis();
                    if (current - lastMsg > 5000) {
                        log.info("{} features processed in {}.", featArray.size(), genome);
                        lastMsg = current;
                    }
                }
                // Write the feature stream.
                try (PrintWriter outStream = new PrintWriter(outFile)) {
                    log.info("Writing feature data to {}.", outFile);
                    outStream.write(featArray.toJson());
                }
                // Write the subsystem stream.
                try (PrintWriter subStream = new PrintWriter(subFile)) {
                    log.info("Writing subsystem data to {}.", subFile);
                    subStream.write(subArray.toJson());
                }
            }
        }
        log.info("{} genomes processed, {} files copied, {} substitutions, {} subsystem records output.",
                gCount, copies, substitutions, subRecords);
    }

    /**
     * Compute the role of a feature in a subsystem.
     *
     * @param sub		subsystem row
     * @param product	function string for feature
     *
     * @return the role for the subsystem
     */
    private String computeRole(SubsystemRow sub, String function) {
        String retVal = null;
        String[] roles = Feature.rolesOfFunction(function);
        List<SubsystemRow.Role> subRoles = sub.getRoles();
        for (int i = 0; i < roles.length && retVal == null; i++) {
            final String name = roles[i];
            boolean found = subRoles.stream().anyMatch(x -> x.getName().equals(name));
            if (found)
                retVal = name;
        }
        return retVal;
    }

}
