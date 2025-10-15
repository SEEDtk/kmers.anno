/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.Option;
import org.theseed.basic.BaseProcessor;
import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.LineReader;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.Role;
import org.theseed.proteins.RoleMap;
import org.theseed.proteins.kmers.KmerReference;
import org.theseed.proteins.kmers.RoleCounter;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.FastaOutputStream;
import org.theseed.sequence.ProteinKmers;
import org.theseed.sequence.Sequence;

/**
 * This method builds a discriminating-kmer database for a fixed set of roles.  It takes as input a role map that can be used to convert role
 * names to role IDs and a file containing a list of the interesting role IDs.  The genomes in a GTO directory are then processed.  The processing
 * is done in two passes.  In the first pass, a kmer database is built from proteins for the interesting roles.  The database is then pruned until
 * only kmers unique to a single role are left.  Proteins that have uninteresting roles are cached in a FASTA file.  The FASTA file is then read
 * back in and used to find kmers that are unique to the interesting roles.  The result is a database that can be used to find interesting roles
 * in genomes that are annotated with a previous, incompatible role-naming system.  The intent is to be able to use this database to create
 * testing sets for the evaluator.
 *
 * The positional parameters are the name of the role definition file, the name of the interesting-role file, and the name of the input GTO
 * directory.
 *
 * The command-line options are as follows:
 *
 * -h	display command-line usage
 * -v	display more detailed progress messages
 * -g	if specified, the name of a file containing genome IDs in the first column; only the identified genomes will be processed
 * -t	the name of a temporary work directory; the default is the current directory
 * -K	protein kmer size
 *
 * @author Bruce Parrello
 *
 */
public class BuildKmerProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(BuildKmerProcessor.class);
    /** map of kmers to role counts */
    private Map<String, RoleCounter> kmerMap;
    /** temporary file for protein FASTA */
    private File fastaBufferFile;
    /** role definition map */
    private RoleMap roleMap;
    /** set of IDs for useful roles */
    private Set<String> goodRoles;
    /** set of IDs for genomes of interest */
    private Set<String> goodGenomes;

    // COMMAND-LINE OPTIONS

    /** optional file of genome IDs to use */
    @Option(name = "-g", aliases = { "--genomes" }, metaVar = "genomeFile.tbl", usage = "file of acceptable genome IDs")
    private File genomeFile;

    /** name of directory for temporary files */
    @Option(name = "-t", aliases = { "--workDir" }, metaVar = "Temp", usage = "temporary file work directory")
    private File workDir;

    /** kmer length */
    @Option(name = "-K", aliases = { "--kmer" }, metaVar = "10", usage = "protein kmer length (default 9)")
    protected void setKmers(int newSize) {
        KmerReference.setKmerSize(newSize);
    }

    /** role definition file */
    @Argument(index = 0, metaVar = "roles.in.subsystems", usage = "role definition file")
    private File roleMapFile;

    /** interesting-role file */
    @Argument(index = 1, metaVar = "roles.to.use", usage = "interesting role file")
    private File roleIdFile;

    /** input GTO directory */
    @Argument(index = 2, metaVar = "genomeDir", usage = "input genome directory")
    private File gtoDir;

    @Override
    protected void setDefaults() {
        this.genomeFile = null;
        this.workDir = new File(System.getProperty("user.dir"));
        this.goodGenomes = null;
        KmerReference.setKmerSize(8);
    }

    @Override
    protected void validateParms() throws IOException {
        // Verify the good-genome list.
        if (this.genomeFile == null)
            log.info("No genome filtering.");
        else {
            if (! this.genomeFile.canRead())
                throw new FileNotFoundException("Good-genome file " + this.genomeFile + " not found or unreadable.");
            this.goodGenomes = TabbedLineReader.readSet(this.genomeFile, "1");
            log.info("{} genome IDs read from genome-filter file.", this.goodGenomes.size());
        }
        // Read in the role map.
        log.info("Reading role definitions from {}.", this.roleMapFile);
        this.roleMap = RoleMap.load(this.roleMapFile);
        // Verify the useful-role list.
        if (! this.roleIdFile.canRead())
            throw new FileNotFoundException("Good-role file " + this.roleIdFile + " not found or unreadable.");
        this.goodRoles = LineReader.readSet(this.roleIdFile);
        // Create the temporary buffer FASTA file.
        if (! this.workDir.isDirectory())
            throw new FileNotFoundException("Temporary-file directory " + this.workDir + " not found or invalid.");
        this.fastaBufferFile = File.createTempFile("buffer", ".faa", this.workDir);
        this.fastaBufferFile.deleteOnExit();
        // Verify the genome directory.
        if (! this.gtoDir.isDirectory())
            throw new FileNotFoundException("Genome directory " + this.gtoDir + " not found or invalid.");    
    }

    @Override
    protected void runCommand() throws Exception {
        // Create the kmer map.
        this.kmerMap = new HashMap<>(this.goodRoles.size() * 700000);
        // We use this to count the number of buffered proteins.
        int bufferedCount = 0;
        // Open the fasta buffer file for output.
        try (FastaOutputStream fastaBuffer = new FastaOutputStream(this.fastaBufferFile)) {
            // Loop through the genomes.
            GenomeDirectory genomes = new GenomeDirectory(this.gtoDir);
            log.info("{} genomes to scan.", genomes.size());
            for (Genome genome : genomes) {
                if (this.goodGenomes == null || this.goodGenomes.contains(genome.getId())) {
                    log.info("Processing {}.", genome);
                    int gBuffCount = 0;
                    int roleCount = 0;
                    int kmerCount = 0;
                    // Run through the pegs.  Good roles are processed for kmers, bad roles are
                    // saved.  A role is only good if it is the sole role of the function, an
                    // artifact of this pipeline's special purpose.
                    for (Feature peg : genome.getPegs()) {
                        List<Role> pegRoles = peg.getUsefulRoles(this.roleMap).stream().filter(x -> this.goodRoles.contains(x.getId())).collect(Collectors.toList());
                        if (pegRoles.isEmpty()) {
                            Sequence pegSequence = new Sequence(peg.getId(), "",
                                    peg.getProteinTranslation());
                            fastaBuffer.write(pegSequence);
                            gBuffCount++;
                            bufferedCount++;
                        } else if (pegRoles.size() == 1) {
                            // Here we have an interesting role.
                            ProteinKmers pegKmers = new ProteinKmers(peg.getProteinTranslation());
                            for (String kmer : pegKmers) {
                                String pegRole = pegRoles.get(0).getId();
                                RoleCounter counter = this.kmerMap.computeIfAbsent(kmer, k -> new RoleCounter(pegRole));
                                boolean good = counter.count(pegRole);
                                if (good) kmerCount++;
                            }
                            roleCount++;
                        }
                    }
                    log.info("{} interesting pegs found, {} buffered, {} good kmers found.", roleCount, gBuffCount, kmerCount);
                }
            }
        }
        // Now we run through the map, deleting bad kmers.
        int deleteCount = 0;
        Iterator<Map.Entry<String, RoleCounter>> iter = this.kmerMap.entrySet().iterator();
        while (iter.hasNext()) {
            Map.Entry<String, RoleCounter> entry = iter.next();
            if (! entry.getValue().isGood()) {
                iter.remove();
                deleteCount++;
            }
        }
        log.info("{} non-unique kmers deleted.  {} discriminating kmers left.  {} proteins buffered.",
                deleteCount, this.kmerMap.size(), bufferedCount);
        // Read back the buffered proteins.  We will count the reads and the kmers deleted.
        bufferedCount = 0;
        deleteCount = 0;
        try (FastaInputStream fastaBuffer = new FastaInputStream(this.fastaBufferFile)) {
            // The actual loop is trivial.  Any kmer we find is deleted.
            for (Sequence seq : fastaBuffer) {
                ProteinKmers kmers = new ProteinKmers(seq.getSequence());
                for (String kmer : kmers) {
                    RoleCounter counter = this.kmerMap.remove(kmer);
                    if (counter != null) deleteCount++;
                }
                bufferedCount++;
                if (log.isInfoEnabled() && bufferedCount % 10000 == 0)
                    log.info("{} proteins processed, {} kmers deleted.", bufferedCount, deleteCount);
            }
        }
        log.info("{} discriminating kmers remaining.", this.kmerMap.size());
        // This will count the kmers found for each role.
        CountMap<String> kmerCounts = new CountMap<>();
        for (Map.Entry<String, RoleCounter> entry : this.kmerMap.entrySet()) {
            String roleId = entry.getValue().getRoleId();
            kmerCounts.count(roleId);
            System.out.println(entry.getKey() + "\t" + roleId);
        }
        if (log.isWarnEnabled()) {
            for (String roleId : this.goodRoles) {
                if (kmerCounts.getCount(roleId) == 0)
                    log.warn("No kmers found for {}: {}.", roleId, this.roleMap.getName(roleId));
            }
        }
    }

}
