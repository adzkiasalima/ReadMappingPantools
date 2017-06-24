/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignerpantools;

import java.io.InputStreamReader;
import org.neo4j.graphdb.GraphDatabaseService;
import genome.SequenceDatabase;
import index.IndexDatabase;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryUsage;
import org.neo4j.graphdb.DynamicLabel;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.RelationshipType;

/**
 *
 * @author adzki001
 */
public class AlignerPantools {

    /**
     * @param args the command line arguments
     */
    
    public static String GRAPH_DATABASE_PATH = "/databases/graph.db/";
    public static String INDEX_DATABASE_PATH = "/databases/index.db/";
    public static String GENOME_DATABASE_PATH = "/databases/genome.db/";
    
    public static GraphDatabaseService graphDb;
    public static IndexDatabase indexDb;
    public static SequenceDatabase genomeDb;
    public static SequenceDatabase sequenceDb;
    public static int MAX_TRANSACTION_SIZE = 100;    //   The number of transactions to be committed in batch
    public static int cores = Runtime.getRuntime().availableProcessors() / 2 + 1;
    
    public static Label pangenome_label = DynamicLabel.label("pangenome");
    public static Label genome_label = DynamicLabel.label("genome");
    public static Label sequence_label = DynamicLabel.label("sequence");
    public static Label node_label = DynamicLabel.label("node");
    public static Label degenerate_label = DynamicLabel.label("degenerate");
//    public static Label panproteome_label = DynamicLabel.label("panproteome");
//    public static Label proteome_label = DynamicLabel.label("proteome");
    public static Label annotation_label = DynamicLabel.label("annotation");
    public static Label variation_label = DynamicLabel.label("variation");
    public static Label gene_label = DynamicLabel.label("gene");
    public static Label coding_gene_label = DynamicLabel.label("coding_gene");
    public static Label RNA_label = DynamicLabel.label("RNA");
    public static Label mRNA_label = DynamicLabel.label("mRNA");
    public static Label tRNA_label = DynamicLabel.label("tRNA");
    public static Label rRNA_label = DynamicLabel.label("rRNA");
    public static Label CDS_label = DynamicLabel.label("CDS");
    public static Label exon_label = DynamicLabel.label("exon");
    public static Label intron_label = DynamicLabel.label("intron");
    public static Label feature_label = DynamicLabel.label("feature");
    public static Label broken_protein_label = DynamicLabel.label("broken_protein");
    public static Label orthology_group_lable = DynamicLabel.label("orthology_group");
    public static Label homology_group_lable = DynamicLabel.label("homology_group");
    public static Label kmer_lable = DynamicLabel.label("kmer");

    
    public static enum RelTypes implements RelationshipType {
        FF, FR, RF, RR,
        has, // for pointing to genome and sequence nodes
        starts,
        stops,
        has_homolog, // for pointing to gene nodes from the homology group
        has_ortholog, // for pointing to gene nodes from the orthology group
        splits_into,
        codes_for,// for connecting genes to mRNAs
        is_parent_of,
        contributes_to,// for connecting CDSs to mRNA
        branches, //to connect tree nodes
        visits,
        precedes, 
        is_homolog_to,
        best_hits
    }

    public static long startTime;
    public static long phaseTime;
    public static int num_nodes;
    public static int num_degenerates;
    public static int num_edges;
    public static int num_bases;
    
    public static void main(String[] args) {
        // TODO code application logic here
        if (args.length < 2 || args[1].equals("--help") || args[1].equals("-h")) {
            print_help_comment();
            System.exit(1);
        }
        
        
    }
    
    private static void print_help_comment() {
        System.out.println("Usage: java - jar "
                + "PATH_TO_THE_JAR_FILE/AlignerPantools.jar "
                + "PATH_TO_THE_PANGENOME_DB "
                + "PATH_TO_THE_READS_FILE" 
             
        );
        
        
    }
    
    
    /**
     * Estimates and prints the peak memory used during the execution of the program. 
     */
    public static void print_peak_memory() {
        long memoryUsage = 0;
        try {
            for (MemoryPoolMXBean pool : ManagementFactory.getMemoryPoolMXBeans()) {
                MemoryUsage peak = pool.getPeakUsage();
                memoryUsage += peak.getUsed();
            }
            System.out.println("Peak memory : " + memoryUsage / 1024 / 1024 + " MB");
        } catch (Throwable t) {
            System.err.println("Exception in agent: " + t);
        }
    }
    
    /**
     * Writes a sequence in a FASTA file with specified length for lines.
     * 
     * @param fasta_file The FASTA file object
     * @param seq The sequence to be written in the file.
     * @param length Length of the lines.
     */    
    public static void write_fasta(BufferedWriter fasta_file, StringBuilder seq, int length) {
        int i;
        try {
            for (i = 1; i <= seq.length(); ++i) {
                fasta_file.write(seq.charAt(i - 1));
                if (i % length == 0) {
                    fasta_file.write("\n");
                }
            }
            fasta_file.write("\n");
        } catch (IOException ioe) {

        }

    }    
    
    /**
     * Return reverse complement of the given string.
     * 
     * @param s The input string
     * @return The reverse complement of the input string
     */     
    public static void reverse_complement(StringBuilder s) {
        char ch;
        int i, j;
        for ( i=0, j = s.length() - 1; i < j; ++i, --j) {
            ch = s.charAt(i);
            s.setCharAt(i, complement(s.charAt(j)));
            s.setCharAt(j, complement(ch));
        }
        if (i == j)
            s.setCharAt(i, complement(s.charAt(i)));
    }

    public static char complement(char ch) {
        switch (ch) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
            case 'R':
                return 'Y';
            case 'Y':
                return 'R';
            case 'K':
                return 'M';
            case 'M':
                return 'K';
            case 'B':
                return 'V';
            case 'V':
                return 'B';
            case 'D':
                return 'H';
            case 'H':
                return 'D';
            default:
                return ch;
        }
    }
    
        /**
     * Executes a shell command. 
     * @param command The command
     * @return The output of the bash command
     */
    public static String executeCommand(String command) {
        StringBuilder exe_output = new StringBuilder();
        String line = "";
        Process p;
        try {
            p = Runtime.getRuntime().exec(command);
            p.waitFor();
            BufferedReader reader = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((line = reader.readLine()) != null) {
                exe_output.append(line + "\n");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return exe_output.toString();
    }    
    
}
