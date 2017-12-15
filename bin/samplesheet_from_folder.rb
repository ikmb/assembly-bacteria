#!/bin/env ruby
# == NAME
# samplesheet_from_folder.rb
#
# == USAGE
# ./this_script.rb [ -h | --help ]
#[ -i | --infile ] |[ -o | --outfile ] | 
# == DESCRIPTION
# A script to produce a basic sample sheet for exome processing from a folder of fastQ files
#
# == OPTIONS
# -h,--help Show help
# -i,--infile=INFILE input file
# -o,--outfile=OUTFILE : output file

#
# == EXPERT OPTIONS
#
# == AUTHOR
#  Marc Hoeppner, mphoeppner@gmail.com

require 'optparse'
require 'ostruct'


### Define modules and classes here


### Get the script arguments and open relevant files
options = OpenStruct.new()
opts = OptionParser.new()
opts.banner = "Reads Fastq files from a folder and writes a sample sheet to STDOUT"
opts.separator ""
opts.on("-f","--folder", "=FOLDER","Folder to scan") {|argument| options.folder = argument }
opts.on("-h","--help","Display the usage information") {
 puts opts
 exit
}

opts.parse! 

abort "Folder not found (#{options.folder})" unless File.directory?(options.folder)

fastq_files = Dir["#{options.folder}/*_R*.fastq.gz"]

groups = fastq_files.group_by{|f| f.split("/")[-1].split(/_R[1,2]/)[0] }

puts "SampleID;libraryID;Date;R1;R2"

#G00076-L2_S19_L003_R1_001.fastq.gz

groups.each do |group, files|

        left,right = files.sort.collect{|f| File.absolute_path(f)}

        library = group.split("_")[0]
        sample = group.split("_")[0]


        puts "#{sample};#{library};#{left};#{right}"

end


