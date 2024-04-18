using FASTX
using ArgParse
using Printf
using Random

include("restriction_enzymes.jl")

struct FragmentInfo
    genome_size::Int
    bases_covered::Int
    num_fragments::Int
    total_fragments_before_filtering::Int
    fragments::Vector{RestrictionFragment}
end

function parse_commandline()

    argument_parser = ArgParseSettings(
        description="Program performs simulated restriction digest.",
        epilog="Erick Samera",
        version="", add_version=true
    )

    @add_arg_table argument_parser begin
        "use-ref"
            help = "perform catalysis on given genome"
            action = :command
        "no-ref"
            help = "perform catalysis on simulated genome of given length and GC content"
            action = :command
    end

    @add_arg_table argument_parser["use-ref"] begin
    "--min-size", "-m"
        help = "minimize size of fragments (bp)"
        arg_type = Int
        default = 300
    "--max-size", "-M"
        help = "maximum size of fragments (bp)"
        arg_type = Int
        default = 600

    "--pretty"
        help = "print prettier output of restriction digest statistics"
        arg_type = Bool
        default = true
    "--csv"
        help = "print csv-formatted output of restriction digest statistics"
        arg_type = Bool
        default = true

    "--gff3"
        help = "print gff3-formatted output of restriction digest fragments"
        arg_type = Bool
        default = false

    "genome_path"
        help = "path to reference genome (.fna/.fasta/.fa/etc.)"
        required = true
    "enzyme_1"
        help = "restriction enzyme 1"
        required = true
    "enzyme_2"
        help = "restriction enzyme 2"
        required = true
    end

    @add_arg_table argument_parser["no-ref"] begin
        "--min-size", "-m"
            help = "minimize size of fragments (bp) to simulate size filtering"
                arg_type = Int
            default = 300
        "--max-size", "-M"
            help = "maximum size of fragments (bp)  to simulate size filtering"
                arg_type = Int
            default = 600

        "genome_size"
            help = "genome size (bp)"
                arg_type = Int
            required = true
        "gc_content"
            help = "GC content (%)"
                arg_type = Float32
            required = true
        "enzyme_1"
            help = "restriction enzyme 1"
            required = true
        "enzyme_2"
            help = "restriction enzyme 2"
            required = true

    "--pretty"
        help = "print prettier output of restriction digest statistics"
        arg_type = Bool
        default = true
    "--csv"
        help = "print csv-formatted output of restriction digest statistics"
        arg_type = Bool
        default = true

    "--seed"
        help = "seed for random sequence generation (-1 = None)"
        arg_type = Int
        default = -1
    end

    return parse_args(argument_parser)
end

function digest_sequence(sequence::String, enzyme1::Enzyme, enzyme2::Enzyme, min_size::Int, max_size::Int, chromosome_id::String)
    initial_fragments = [RestrictionFragment(sequence, "original", "original", 1, length(sequence), chromosome_id)]
    
    fragments_after_first_digest = catalyze(enzyme1, initial_fragments)

    bases_covered = 0
    total_fragments = 0
    num_fragments_after_filter = 0

    fragments_collected = []

    for fragment in fragments_after_first_digest
        fragments_after_second_digest = catalyze(enzyme2, [fragment])
        
        for fragment2 in fragments_after_second_digest
            if fragment2.left_end != fragment2.right_end
                len = length(fragment2.sequence)
                total_fragments += 1
                if min_size < len < max_size
                    bases_covered += len
                    num_fragments_after_filter += 1
                    push!(fragments_collected, fragment2)
                end
            end
        end
    end
    
    return bases_covered, num_fragments_after_filter, total_fragments, fragments_collected
end

function print_csv(info::FragmentInfo, parsed_args::Dict{String, Any})
    @printf("%s-%s,%s,%s,%s,%s\n", parsed_args[parsed_args["%COMMAND%"]]["enzyme_1"], parsed_args[parsed_args["%COMMAND%"]]["enzyme_2"], info.bases_covered, info.num_fragments, info.total_fragments_before_filtering, info.bases_covered / info.genome_size * 100)
end

function print_pretty(info::FragmentInfo, parsed_args::Dict{String, Any})
    @printf("# %s-%s\n", parsed_args[parsed_args["%COMMAND%"]]["enzyme_1"], parsed_args[parsed_args["%COMMAND%"]]["enzyme_2"])
    println("# Genome size: ", info.genome_size)
    println("# Bases covered: ", info.bases_covered)
    println("# Number of fragments after filtering: ", info.num_fragments)
    println("# Number of fragments before filtering: ", info.total_fragments_before_filtering)
    println("# Percentage of bases covered: ", info.bases_covered / info.genome_size * 100, "%")
end

function generate_sequence(length::Int, gc_content::Float32)
    num_gc = Int(round(length * (gc_content / 100)))
    num_at = length - num_gc

    sequence = vcat(fill('G', num_gc ÷ 2), 
                    fill('C', num_gc ÷ 2), 
                    fill('A', num_at ÷ 2), 
                    fill('T', num_at ÷ 2))

    if num_gc % 2 != 0
        push!(sequence, rand(['G', 'C']))
    end

    if num_at % 2 != 0
        push!(sequence, rand(['A', 'T']))
    end
    return join(shuffle!(sequence))
end

function main_ref_digestion(parsed_args)

    collected_fragments = []

    open(FASTA.Reader, parsed_args[parsed_args["%COMMAND%"]]["genome_path"]) do reader
        total_genome_size = 0
        total_bases_covered = 0
        total_num_fragments = 0
        total_fragments_before_filtering = 0
        
        collected_fragments = []
        
        for record in reader
            sequence = String(FASTX.sequence(record))
            total_genome_size += length(sequence)
    
            bases_covered, num_fragments, total_fragments, fragments = digest_sequence(
                sequence, 
                enzyme_dict[parsed_args[parsed_args["%COMMAND%"]]["enzyme_1"]], 
                enzyme_dict[parsed_args[parsed_args["%COMMAND%"]]["enzyme_2"]], 
                parsed_args[parsed_args["%COMMAND%"]]["min-size"],
                parsed_args[parsed_args["%COMMAND%"]]["max-size"],
                String(FASTX.identifier(record))
            )

            total_bases_covered += bases_covered
            total_num_fragments += num_fragments
            total_fragments_before_filtering += total_fragments
            append!(collected_fragments, fragments)
        end
        
        info = FragmentInfo(total_genome_size, total_bases_covered, total_num_fragments, total_fragments_before_filtering, collected_fragments)
        
        if parsed_args[parsed_args["%COMMAND%"]]["gff3"]
            println("##gff-version 3")
            for fragment in info.fragments
                @printf("%s\t%s\t%s\t%d\t%d\t.\t.\t.\tID=fragment_%d;Name=fragment_%d\n", 
                        fragment.id, # Using the chromosome or contig ID from the fragment
                        "simRAD", 
                        "restriction_fragment",
                        fragment.start_pos,
                        fragment.end_pos,
                        fragment.start_pos, # Just using start_pos as a unique ID here
                        fragment.start_pos) 
            end
        end

        if parsed_args[parsed_args["%COMMAND%"]]["csv"]
            print_csv(info, parsed_args)
        end
        if parsed_args[parsed_args["%COMMAND%"]]["pretty"]
            print_pretty(info, parsed_args)
        end
    end
end

function main_non_ref_digestion(parsed_args)

    if parsed_args[parsed_args["%COMMAND%"]]["seed"] != -1
        Random.seed!(parsed_args[parsed_args["%COMMAND%"]]["seed"])
    end

    sequence = generate_sequence(parsed_args[parsed_args["%COMMAND%"]]["genome_size"], parsed_args[parsed_args["%COMMAND%"]]["gc_content"])

    total_genome_size = 0
    total_bases_covered = 0
    total_num_fragments = 0
    total_fragments_before_filtering = 0
    
    total_genome_size += length(sequence)

    bases_covered, num_fragments, total_fragments = digest_sequence(
        sequence, 
        enzyme_dict[parsed_args[parsed_args["%COMMAND%"]]["enzyme_1"]], 
        enzyme_dict[parsed_args[parsed_args["%COMMAND%"]]["enzyme_2"]], 
        parsed_args[parsed_args["%COMMAND%"]]["min-size"], 
        parsed_args[parsed_args["%COMMAND%"]]["max-size"],
        "simulated_genomer"
    )

    total_bases_covered += bases_covered
    total_num_fragments += num_fragments
    total_fragments_before_filtering += total_fragments
    
    info = FragmentInfo(total_genome_size, total_bases_covered, total_num_fragments, total_fragments_before_filtering, [])

    if parsed_args[parsed_args["%COMMAND%"]]["csv"]
        print_csv(info, parsed_args)
    end
    if parsed_args[parsed_args["%COMMAND%"]]["pretty"]
        print_pretty(info, parsed_args)
    end
end

function main()

    parsed_args = parse_commandline()

    if !haskey(enzyme_dict, parsed_args[parsed_args["%COMMAND%"]]["enzyme_1"])
        println("Enzyme 1 is not recognized!")
    end

    if !haskey(enzyme_dict, parsed_args[parsed_args["%COMMAND%"]]["enzyme_2"])
        println("Enzyme 2 is not recognized!")
    end

    println(stderr, "## ", parsed_args[parsed_args["%COMMAND%"]]["enzyme_1"], "-",parsed_args[parsed_args["%COMMAND%"]]["enzyme_2"])

    if parsed_args["%COMMAND%"] == "use-ref"
        main_ref_digestion(parsed_args)
    end

    if parsed_args["%COMMAND%"] == "no-ref"
        main_non_ref_digestion(parsed_args)
    end
end

main()