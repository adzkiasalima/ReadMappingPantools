/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;

import java.io.File;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import static alignerpantools.AlignerPantools.GRAPH_DATABASE_PATH;
import static alignerpantools.AlignerPantools.INDEX_DATABASE_PATH;
import static alignerpantools.AlignerPantools.GENOME_DATABASE_PATH;
import static alignerpantools.AlignerPantools.MAX_TRANSACTION_SIZE;
import alignerpantools.AlignerPantools.RelTypes;
import static alignerpantools.AlignerPantools.complement;
import static alignerpantools.AlignerPantools.degenerate_label;
import static alignerpantools.AlignerPantools.genomeDb;
import static alignerpantools.AlignerPantools.graphDb;
import static alignerpantools.AlignerPantools.indexDb;
import static alignerpantools.AlignerPantools.node_label;
import static alignerpantools.AlignerPantools.num_bases;
import static alignerpantools.AlignerPantools.num_degenerates;
import static alignerpantools.AlignerPantools.num_edges;
import static alignerpantools.AlignerPantools.num_nodes;
import static alignerpantools.AlignerPantools.pangenome_label;
import static alignerpantools.AlignerPantools.sequence_label;
import static alignerpantools.AlignerPantools.startTime;
import genome.SequenceDatabase;
import index.IndexDatabase;
import index.IndexPointer;
import index.kmer;
import java.util.Arrays;
import org.neo4j.graphdb.Direction;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;

/**
 *
 * @author adzki001
 */
public class readsLayer {
    /**
     * Adds new genomes to an available pangenome.
     * 
     * @param genome_paths_file Path to the FASTA genome files. 
     * @param pangenome_path Path to the database folder
     */
    
    private kmer fwd_kmer, rev_kmer, k_mer;
    private long seq_len;
    private static int K, genome, sequence, position;
    private long curr_index;
    private byte curr_side;
    private Node curr_node;
    private Node db_node;
    private Node new_node;
    private Node degenerate_node;
    private IndexPointer pointer;
    private int fwd_code, rev_code;
    private boolean finish;
    private static int ANCHOR_DISTANCES = 100; // Distance between two anchor nodes
    
    public void add_genomes(String genome_paths_file, String pangenome_path) {
        int j, len, previous_num_genomes;
        long byte_number = 0;
        // address : iterator 
        int[] address = new int[4];  
        // start : 
        // seq_node : 
        Node start, seq_node;
        if (new File(pangenome_path + GRAPH_DATABASE_PATH).exists()) {
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(pangenome_path + GRAPH_DATABASE_PATH))
                    .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
            registerShutdownHook(graphDb);
            startTime = System.currentTimeMillis();
            try (Transaction tx = graphDb.beginTx()) {
                db_node = graphDb.findNodes(pangenome_label).next();
                if (db_node == null) {
                    System.out.println("Can not locate database node!");
                    System.exit(1);
                }
            // Reads the properties of the pangenome    
                K = (int) db_node.getProperty("k_mer_size");
                num_nodes = (int) db_node.getProperty("num_nodes");
                num_edges = (int) db_node.getProperty("num_edges");
                num_degenerates = (int) db_node.getProperty("num_degenerate_nodes");
                num_bases = (int) db_node.getProperty("num_bases");
                previous_num_genomes = (int) db_node.getProperty("num_genomes");
            // if the genome database is not available, reconstruct it.    
                if (!Files.exists(Paths.get(pangenome_path + GENOME_DATABASE_PATH))) 
                {
                // read genomes information from the graph and rebuild the genomes database
                    genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH, graphDb);
                    // create a new instance seq
                    StringBuilder seq = new StringBuilder();
                    for (address[0] = 1; address[0] <= genomeDb.num_genomes; ++address[0]) {
                        for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                            seq_node = graphDb.findNode(sequence_label, "number", address[0] + "_" + address[1]);
                            start = seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                            address[2] = 1;
                            address[3] = (int) genomeDb.sequence_length[address[0]][address[1]];
                            extract_sequence(seq, new IndexPointer(start.getId(), true, 0, -1l), address);
                            len = seq.length();
                            if (len % 2 == 1) {
                                --len;
                            }
                            for (j = 0; j < len; j += 2, ++byte_number) {
                                genomeDb.genomes_buff[(int) (byte_number / genomeDb.parts_size[0])].put((byte) ((genomeDb.binary[seq.charAt(j)] << 4) | genomeDb.binary[seq.charAt(j + 1)]));
                            }
                            if (len == seq.length() - 1) {
                                genomeDb.genomes_buff[(int) (byte_number / genomeDb.parts_size[0])].put((byte) (genomeDb.binary[seq.charAt(len)] << 4));
                                ++byte_number;
                            }
                        }
                    }
                    genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH); //Readable only
                } else {
                    genomeDb = new SequenceDatabase(pangenome_path + GENOME_DATABASE_PATH);
                }
                genomeDb.add_genomes(pangenome_path + GENOME_DATABASE_PATH, genome_paths_file);
                indexDb = new IndexDatabase(pangenome_path + INDEX_DATABASE_PATH, genome_paths_file, genomeDb, graphDb, previous_num_genomes);
                tx.success();
            }
        // the sequences should be dropped out as they will change and add_sequence_properties() function will rebuild them.    
            drop_nodes_property("sequence");
        // the edge colors should be dropped out as they will change and localize_nodes() function will rebuild them again.    
            drop_edges_colors();
            construct_pangenome(previous_num_genomes);
            System.out.println("Number of kmers:   " + indexDb.length());
            System.out.println("Number of nodes:   " + num_nodes);
            System.out.println("Number of edges:   " + num_edges);
            System.out.println("Number of bases:   " + num_bases);
            System.out.println("Number of degenerate nodes:   " + num_degenerates);
            try (Transaction tx = graphDb.beginTx()) {
                db_node.setProperty("k_mer_size", K);
                db_node.setProperty("num_k_mers", indexDb.length());
                db_node.setProperty("num_nodes", num_nodes);
                db_node.setProperty("num_degenerate_nodes", num_degenerates);
                db_node.setProperty("num_edges", num_edges);
                db_node.setProperty("num_genomes", genomeDb.num_genomes);
                db_node.setProperty("num_bases", num_bases);
                tx.success();
            }
            graphDb.shutdown();
            genomeDb.close();
            indexDb.close();
            File directory = new File(pangenome_path + GRAPH_DATABASE_PATH);
            for (File f : directory.listFiles()) {
                if (f.getName().startsWith("neostore.transaction.db.")) {
                    f.delete();
                }
            }
            System.out.println("graph.db size: " + getFolderSize(new File(pangenome_path + GRAPH_DATABASE_PATH)) + " MB");
            System.out.println("index.db size: " + getFolderSize(new File(pangenome_path + INDEX_DATABASE_PATH)) + " MB");
            System.out.println("genome.db size: " + getFolderSize(new File(pangenome_path + GENOME_DATABASE_PATH)) + " MB");
        } else {
            System.out.println("No database found in " + pangenome_path);
            System.exit(1);
        }
    }
    
    /**
     * Shuts down the graph database if the program halts unexpectedly.
     * 
     * @param graphDb The graph database object 
     */
    private static void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }

    /**
     * Extracts the genomic region belonging to the specified sequence starting at th specified node.
     * 
     * @param seq Will contains the sequence after function ends.
     * @param start_ptr A pangenome pointer which points to the node where the sequence starts.
     * @param address An array determining {genome, sequence, begin, end}properties of the sequence.
     */
    public static void extract_sequence(StringBuilder seq, IndexPointer start_ptr, int[] address) {
        Relationship rel;
        Node neighbor, node;
        int[] addr = new int[]{address[0],address[1],address[2],address[3]};
        int begin = addr[2] - 1, end = addr[3] - 1;
        int loc, node_len, neighbor_len, seq_len, position;
        String rel_name;
        seq_len = end - begin + 1;
        seq.setLength(0);
        loc = begin;
        position=start_ptr.offset;
        node=graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
    // Takes the part of the region lies in the first node of the path that region takes in the graph    
        if (start_ptr.canonical) {
            if (position + seq_len - 1 <= node_len - 1) { // The whole sequence lies in this node
                loc += append_fwd(seq, (String) node.getProperty("sequence"), position, position + seq_len - 1);
            } else {
                loc += append_fwd(seq, (String) node.getProperty("sequence"), position, node_len - 1);
            }
        } else {
            if (position - (seq_len - 1) >= 0) { // The whole sequence lies in this node
                loc += append_rev(seq, (String) node.getProperty("sequence"), position - (seq_len - 1), position);
            } else {
                loc += append_rev(seq, (String) node.getProperty("sequence"), 0, position);
            }
        }
    //  traverse the path of the region   
        while (seq.length() < seq_len ) {
            addr[2] = loc - K + 1;
            rel = get_outgoing_edge(node, addr);
            neighbor = rel.getEndNode();
            rel_name = rel.getType().name();
            neighbor_len = (int) neighbor.getProperty("length");
            if (rel_name.charAt(1) == 'F') // Enterring forward side
                if (seq.length() + neighbor_len - K + 1 > seq_len) // neighbor is the last node of the path
                    //loc += append_fwd(seq, (String) neighbor.getProperty("sequence"), 0, seq_len - seq.length() + K - 2);
                    loc += append_fwd(seq, (String) neighbor.getProperty("sequence"), K - 1, seq_len - seq.length() + K - 2);
                else 
                    //loc += append_fwd(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - 1);
                    loc += append_fwd(seq, (String) neighbor.getProperty("sequence"), K - 1, neighbor_len - 1);
            else // Enterring reverse side
                if (seq.length() + neighbor_len - K + 1 > seq_len) // neighbor is the last node of the pat
                    //loc += append_rev(seq, (String) neighbor.getProperty("sequence"), neighbor_len - (seq_len - seq.length()), neighbor_len - 1);
                    loc += append_rev(seq, (String) neighbor.getProperty("sequence"), neighbor_len - K - (seq_len - seq.length()) + 1, neighbor_len - K);
                else 
                    //loc += append_rev(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - 1);
                    loc += append_rev(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K);
            node = neighbor;
        } // while
    }
    
    /**
     * Give the next node of the path to be traversed through.
     * 
     * @param current_node The current node of the path we are located at. 
     * @param address An array which determine the genome, sequence and position of the desirable outgoing edge.
     * @return The outgoing edge.
     */
    public static Relationship get_outgoing_edge(Node current_node, int[] address) {
        String origin;
        int[] occurrence;
        int genome = address[0], sequence = address[1];
        origin = "a" + genome + "_" + sequence;
        for (Relationship r_out : current_node.getRelationships(Direction.OUTGOING)) {
            occurrence = (int[])r_out.getProperty(origin, null);
            if (occurrence != null) {
                if (Arrays.binarySearch(occurrence, address[2])>=0)
                    return r_out;
            }
        }
        return null;
    }
    
    /**
     * Appends substring s[from..to] to the string builder.
     * @param seq String builder.
     * @param s The string.
     * @param from Start position of the substring.
     * @param to Stop position of the substring.
     * @return The length of substring appended to the string builder.
     */
    public static int append_fwd(StringBuilder seq, String s, int from, int to) {
        for (int i = from; i <= to; ++i)
            seq.append(s.charAt(i));
        return to - from + 1;
    }
    
    /**
     * Appends the reverse complement of substring s[from..to] to the string builder.
     * @param seq String builder.
     * @param s The string.
     * @param from Start position of the substring.
     * @param to Stop position of the substring.
     * @return The length of substring appended to the string builder.
     */
    public static int append_rev(StringBuilder seq, String s, int from, int to) {
        for (int i = to; i >= from; --i)
            seq.append(complement(s.charAt(i)));
        return to - from + 1;
    }
    
    /**
     * Returns a pangenome pointer pointing to the specified genomic position.
     * 
     * @param address An integer array lile {genome_number, sequence_number, begin_position, end_position}
     * @return A pointer to the genomic position in the pangenome
     */
    public static IndexPointer locate(int[] addr, int K) {
        int node_start_pos, low, high, mid , node_len, genomic_pos;
        boolean forward;
        Node node, neighbor, seq_node;
        Relationship rel;
        String anchor_sides;
        long[] anchor_nodes;
        int[] anchor_positions;
        int[] address = Arrays.copyOf(addr,addr.length);
        genomic_pos = address[2] - 1;
        seq_node = graphDb.findNode(sequence_label, "number", address[0]+"_"+address[1]);
        anchor_nodes = (long[]) seq_node.getProperty("anchor_nodes");
        anchor_positions = (int[]) seq_node.getProperty("anchor_positions");
        anchor_sides = (String) seq_node.getProperty("anchor_sides");
    // Find the immediate preceding anchor_node, searching in the sorted array of anchor positions.      
        for (low = 0, high = anchor_sides.length() - 1, mid = (low + high) / 2; low <= high; mid = (low + high) / 2) {
            if (genomic_pos < anchor_positions[mid]) {
                high = mid - 1;
            } else if (genomic_pos > anchor_positions[mid]) {
                low = mid + 1;
            } else {
                break;
            }
        }
        if (genomic_pos < anchor_positions[mid]) {
            --mid;
        }
        forward = anchor_sides.charAt(mid) == 'F';
        try (Transaction tx = graphDb.beginTx()) {
            node = graphDb.getNodeById(anchor_nodes[mid]);
            node_start_pos = anchor_positions[mid];
            node_len = (int) node.getProperty("length");
        // Traverse the pangenome from the anchor node until reach to the target
            while (node_start_pos + node_len <= genomic_pos) 
            {
                address[2] = node_start_pos + node_len - K + 1;
                rel = get_outgoing_edge(node, address);
                if (rel == null){
                    System.out.println("Failed to locate address : " + address[0] + " " + address[1] + " "+ address[2]);
                    break;
                }
                neighbor = rel.getEndNode();
                forward = rel.getType().name().charAt(1) == 'F';
                node_start_pos += node_len - K + 1;
                node = neighbor;
                node_len = (int) node.getProperty("length");
            }
            tx.success();
        }
        return new IndexPointer(node.getId(), forward, forward ? genomic_pos - node_start_pos : node_len - 1 - (genomic_pos - node_start_pos), -1l);
    }

    /**
     * Remove the given property from all the nodes and degenerates.
     * @param property 
     */    
    void drop_nodes_property(String property) {
        int i;
        byte[] sequence;
        ResourceIterator<Node> nodes;
        Node node;
        try (Transaction tx = graphDb.beginTx()) {
            nodes = graphDb.findNodes(node_label);
            tx.success();
        }
        while (nodes.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && nodes.hasNext(); ++i) {
                    node = nodes.next();
                    node.removeProperty(property);
                }
                tx.success();
            }
        }
        nodes.close();
        try (Transaction tx = graphDb.beginTx()) {
            nodes = graphDb.findNodes(degenerate_label);
            tx.success();
        }
        while (nodes.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && nodes.hasNext(); ++i) {
                    node = nodes.next();
                    node.removeProperty(property);
                }
                tx.success();
            }
        }
        nodes.close();
    }



    /**
     * Remove the occurrence arrays of the edges.
     */    
    void drop_edges_colors() {
        int i;
        ResourceIterator<Relationship> rels;
        Relationship r;
        try (Transaction tx = graphDb.beginTx()) {
            rels = graphDb.getAllRelationships().iterator();
            tx.success();
        }
        while (rels.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && rels.hasNext(); ++i) {
                    r = rels.next();
                    if (r.isType(RelTypes.FF) || r.isType(RelTypes.FR) || r.isType(RelTypes.RF) || r.isType(RelTypes.RR))
                        for(String p:r.getPropertyKeys())
                            r.removeProperty(p);
                }
                tx.success();
            }
        }
        rels.close();
    }

    private void construct_pangenome(int previous_num_genomes) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Returns size of a given folder.
     * 
     * @param dir The folder File object.
     * @return Size of the folder in MegaBytes
     */
    public static long getFolderSize(File dir) {
        long size = 0;
        for (File file : dir.listFiles()) {
            if (file.isFile()) {
                // System.out.println(file.getName() + " " + file.length());
                size += file.length();
            } else {
                size += getFolderSize(file);
            }
        }
        return size / 1048576 + 1;
    }
}
