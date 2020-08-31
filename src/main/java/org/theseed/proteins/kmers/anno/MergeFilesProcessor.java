/**
 *
 */
package org.theseed.proteins.kmers.anno;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.theseed.io.LineReader;
import org.theseed.utils.BaseProcessor;

/**
 * This command is run after the "build" and "apply" commands.  At this point, we have a testing set file and a training set file.
 * The two files must be merged, with the testing set at the top.  All roles that do not occur in the testing set will be removed
 * from the result file and the roles.to.use file.
 *
 * The positional parameter is the name of the evaluation directory.  It must contain the roles.to.use file, the testing.tbl file,
 * and the training.tbl file.  The modified files will first be backed up to the Backup subdirectory.
 *
 * The command-line options are as follows.
 *
 * -h	command-line usage
 * -v	display more frequent progress messages
 *
 *
 * @author Bruce Parrello
 *
 */
public class MergeFilesProcessor extends BaseProcessor {

    // FIELDS
    /** logging facility */
    protected static Logger log = LoggerFactory.getLogger(MergeFilesProcessor.class);
    /** roles-to-use file */
    private File rolesFile;
    /** testing set file */
    private File testingFile;
    /** training set file */
    private File trainingFile;
    /** column flags; TRUE if the column should be kept */
    private boolean[] keepFlags;
    /** buffer for building output lines */
    private StringBuilder buffer;


    // COMMAND-LINE OPTIONS

    /** evaluation directory */
    @Argument(index = 0, metaVar = "evalDir", usage = "evaluation directory")
    private File evalDir;

    @Override
    protected void setDefaults() {
    }

    @Override
    protected boolean validateParms() throws IOException {
        // Verify the evaluation directory.
        if (! this.evalDir.isDirectory())
            throw new FileNotFoundException("Evaluation directory " + this.evalDir + " not found or invalid");
        // Verify the backup directory.
        File backupDir = new File(this.evalDir, "Backup");
        if (! backupDir.isDirectory()) {
            log.info("Creating backup directory " + backupDir + ".");
            FileUtils.forceMkdir(backupDir);
        }
        // Verify the files.
        this.rolesFile = new File(this.evalDir, "roles.to.use");
        this.testingFile = new File(this.evalDir, "testing.tbl");
        this.trainingFile = new File(this.evalDir, "training.tbl");
        if (! this.rolesFile.canRead())
            throw new FileNotFoundException("Roles-to-use file " + this.rolesFile + " not found or unreadable.");
        if (! this.testingFile.canRead())
            throw new FileNotFoundException("Testing file " + this.testingFile + " not found or unreadable.");
        if (! this.trainingFile.canRead())
            throw new FileNotFoundException("Training file " + this.trainingFile + " not found or unreadable.");
        // Backup the ones we're modifying.
        FileUtils.copyFileToDirectory(this.rolesFile, backupDir);
        FileUtils.copyFileToDirectory(this.trainingFile, backupDir);
        return true;
    }

    @Override
    protected void runCommand() throws Exception {
        // First, we get the training records.
        log.info("Processing training file {}.", this.trainingFile);
        List<String[]> trainLines = new ArrayList<String[]>(1000);
        try (LineReader trainingStream = new LineReader(this.trainingFile)) {
            for (String line : trainingStream)
                trainLines.add(StringUtils.split(line, '\t'));
        }
        this.keepFlags = new boolean[trainLines.get(0).length];
        log.info("{} records and {} columns in training set.", trainLines.size(), keepFlags.length);
        // Now process the testing records.
        List<String[]> testLines = new ArrayList<String[]>(200);
        try (LineReader testingStream = new LineReader(this.testingFile)) {
            for (String line : testingStream) {
                String[] fields = StringUtils.split(line, '\t');
                for (int i = 0; i < keepFlags.length; i++)
                    if (! fields[i].contentEquals("0")) keepFlags[i] = true;
                testLines.add(fields);
            }
        }
        log.info("{} lines read from testing set.", testLines.size());
        if (log.isInfoEnabled()) {
            int colsKept = 0;
            for (int i = 0; i < this.keepFlags.length; i++)
                if (this.keepFlags[i]) colsKept++;
            log.info("{} columns will be kept.", colsKept);
        }
        // Create the output string builder.
        this.buffer = new StringBuilder(10 + 5 * this.keepFlags.length);
        // Now we need to rebuild the training file.
        try (PrintStream trainingStream = new PrintStream(this.trainingFile)) {
            // Write the header line.
            this.writeLine(trainingStream, trainLines.get(0));
            // Write the testing lines.
            log.info("Writing testing set.");
            for (String[] line : testLines)
                this.writeLine(trainingStream, line);
            // Write the training lines.
            for (int i = 1; i < trainLines.size(); i++)
                this.writeLine(trainingStream, trainLines.get(i));
        }
        // Now we rebuild roles.to.use.
        List<String> roleLines = new ArrayList<String>(keepFlags.length);
        log.info("Reading role file.");
        try (LineReader roleStream = new LineReader(this.rolesFile)) {
            int i = 1;
            for (String line : roleStream) {
                if (keepFlags[i])
                    roleLines.add(line);
                i++;
            }
        }
        log.info("Updating role file. {} roles will be kept.", roleLines.size());
        try (PrintStream roleStream = new PrintStream(this.rolesFile)) {
            for (String line : roleLines)
                roleStream.println(line);
        }
        log.info("Files created.");
    }

    /**
     * Write a line to the training stream.  Only the columns indicated by the keep flags are output.
     *
     * @param trainingStream	output stream
     * @param strings			array of fields in this line
     */
    private void writeLine(PrintStream trainingStream, String[] strings) {
        this.buffer.setLength(0);
        // We always keep the first column.
        this.buffer.append(strings[0]);
        for (int i = 1; i < this.keepFlags.length; i++) {
            if (this.keepFlags[i])
                this.buffer.append('\t').append(strings[i]);
        }
        trainingStream.println(buffer.toString());
    }

}
