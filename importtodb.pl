#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use File::Spec;
use File::Basename;
use DBI;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use threads;
use Thread::Queue;
use DateTime;
use POSIX qw( ceil );
use lib '/home/modupe/SCRIPTS/SUB';
use routine;
use passw;

# #CREATING LOG FILES
my $std_out = '/home/modupe/.LOG/RavTAD-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupe/.LOG/RavTAD-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';
my $jobid = "RavenTAD-".`date +%m-%d-%y_%T`;

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";

my ($VERSION, $DATE, $AUTHOR, $CONTACT, $EMAIL) = DEFAULTS();

#--------------------------------------------------------------------------------

my ($help, $manual, $sendemail, $folder2import, $mappingtool,$ardea, $ibis);
our ($verbose, $nosql, $vnosql, $gnosql, $cnosql, $log, $transaction);
our ($metadata, $tab, $excel, $datadb, $gene, $variant, $all, $vep, $annovar, $delete); #command options
our ($file2consider,$connect); #connection and file details
my ($sth,$dbh,$schema); #connect to database;

#data2db options
our ($found);
our (@allgeninfo);
our ($str, $ann, $ref, $seq,$allstart, $allend) = (0,0,0,0,0,0); #for log file
our ($refgenome, $refgenomename, $stranded, $sequences, $annotationfile, $mparameters, $gparameters, $cparameters, $vparameters); #for annotation file
my $additional;

#genes import
our ($bamfile, $fastqcfolder, $alignfile, $staralignfile, $version, $readcountfile, $starcountfile, $genesfile, $deletionsfile, $insertionsfile, $transcriptsgtf, $junctionsfile, $variantfile, $vepfile, $annofile);
our ($kallistofile, $kallistologfile, $salmonfile, $salmonlogfile);
our ($totalreads, $mapped, $alignrate, $deletions, $insertions, $junctions, $genes, $annversion, $diffexpress, $counttool);
my (%ARFPKM,%CHFPKM, %BEFPKM, %CFPKM, %DFPKM, %TPM, %cfpkm, %dfpkm, %tpm, %DHFPKM, %DLFPKM, %dhfpkm, %dlfpkm, %ALL);
my (%HASHDBVARIANT, @VAR, @threads, $queue);
#variant import
our ( %VCFhash, %DBSNP, %extra, %VEPhash, %ANNOhash );
our ($varianttool, $verd, $variantclass, $vcount);
our ($itsnp,$itindel,$itvariants) = (0,0,0);

#nosql append
my (@nosqlrow, $showcase);
#date
my $date = `date +%Y-%m-%d`;

#--------------------------------------------------------------------------------

# -----------------------------------
# CREATING EMAIL NOTIFICATION
if ($sendemail) { NOTIFICATION("Starting Job",0,0); }
# -----------------------------------

# CONNECT TO THE DATABASE
processArguments();

if ($folder2import =~ /\w.*_(\d+)/){
  $dbh = mysql();
  my $libraryidnumber = $1;
  print "JOB:\tImporting Transcriptome analysis Information => $libraryidnumber\n";
  `find $folder2import` or pod2usage ("ERROR: Can not locate \"$folder2import\"");
  opendir (DIR, $folder2import) or pod2usage ("ERROR: $folder2import is not a folder, please specify your sample directory location "); close (DIR);
  my @foldercontent = split("\n", `find $folder2import -type f -print0 | xargs -0 ls -tr `); #get details of the folder
		
  foreach (grep /\.gtf/, @foldercontent) { unless (`head -n 3 $_ | wc -l` <= 0 && $_ =~ /skipped/) { $transcriptsgtf = $_; } }
	$fastqcfolder = (grep /fastqc.zip$/, @foldercontent)[0]; unless ($fastqcfolder) { $fastqcfolder = (grep /fastqc_data.txt$/, @foldercontent)[0]; }
	$alignfile = (grep /summary.txt/, @foldercontent)[0];
  $genesfile = (grep /genes.fpkm/, @foldercontent)[0];
  $deletionsfile = (grep /deletions.bed/, @foldercontent)[0];
  $insertionsfile = (grep /insertions.bed/, @foldercontent)[0];
  $junctionsfile = (grep /junctions.bed/, @foldercontent)[0];
	$bamfile = (grep /.bam$/, @foldercontent)[0];
  $variantfile = (grep /.vcf$/, @foldercontent)[0]; 
  $vepfile = (grep /.vep.txt$/, @foldercontent)[0];
  $annofile = (grep /anno.txt$/, @foldercontent)[0];
	$readcountfile = (grep /.counts$/, @foldercontent)[0];
	$starcountfile = (grep /ReadsPerGene.out.tab$/, @foldercontent)[0];
	$kallistofile = (grep /.tsv$/, @foldercontent)[0];
	$kallistologfile = (grep /run_info.json/, @foldercontent)[0];
	$salmonfile = (grep /.sf$/, @foldercontent)[0];
	$salmonlogfile = (grep /cmd_info.json$/, @foldercontent)[0];
	$staralignfile = (grep /Log.final.out$/, @foldercontent)[0];

	$sth = $dbh->prepare("select libraryid from BirdLibraries where libraryid = $libraryidnumber"); $sth->execute(); $found = $sth->fetch();
  if ($found) { # if sample is not in the database    
    $sth = $dbh->prepare("select libraryid from MappingStats where libraryid = $libraryidnumber"); $sth->execute(); $found = $sth->fetch();
    LOGFILE();
		unless ($found) { 
			#open alignment summary file
      unless ($kallistologfile || $salmonlogfile) {
				if ($alignfile) {
					`head -n 1 $alignfile` =~ /^(\d+)\sreads/; $totalreads = $1;
					open(ALIGN,"<", $alignfile) or die "\nFAILED:\t Can not open Alignment summary file '$alignfile'\n";
					while (<ALIGN>){
						chomp;
						if (/Input/){my $line = $_; $line =~ /Input.*:\s+(\d+)$/;$totalreads = $1;}
						if (/overall/) { my $line = $_; $line =~ /^(\d+.\d+)%\s/; $alignrate = $1;}
						if (/overall read mapping rate/) {
							if ($mappingtool){
								unless ($mappingtool =~ /TopHat/i){
									die "\nERROR:\t Inconsistent Directory Structure, $mappingtool SAM file with TopHat align_summary.txt file found\n";
								}
							} else { $mappingtool = "TopHat"; }
						}
						if (/overall alignment rate/) {
							if ($mappingtool){
								unless ($mappingtool =~ /hisat/i){
									die "\nERROR:\t Inconsistent Directory Structure, $mappingtool LOG file with HISAT align_summary.txt file found\n";
								}
							} else { $mappingtool = "HISAT";}
						}
					} close ALIGN;
					$mapped = ceil($totalreads * $alignrate/100);
				} elsif ($staralignfile && $mappingtool =~ /star/i) {
					`grep "Number of input reads" $staralignfile` =~ /\s(\d+)$/; $totalreads = $1;
					`grep "Uniquely mapped reads %" $staralignfile` =~ /\s(\S+)\%$/; $alignrate = $1;
				} else {
					if ($mappingtool =~ /star/) {die "\nFAILED:\t Can not find STAR Alignment summary file as 'Log.final.out'\n";}
					else {die "\nFAILED:\t Can not find Alignment summary file as 'align_summary.txt'\n";}
				}
			}
     				
			$deletions = undef; $insertions = undef; $junctions = undef;
			if ($deletionsfile){ $deletions = `cat $deletionsfile | wc -l`; $deletions--; } 
			if ($insertionsfile){ $insertions = `cat $insertionsfile | wc -l`; $insertions--; }
			if ($junctionsfile){ $junctions = `cat $junctionsfile | wc -l`; $junctions--; }
						
			#INSERT INTO DATABASE:
			#MapStats table
			if ($mappingtool) { print "NOTICE:\t Importing $mappingtool alignment information for $libraryidnumber to MappingStats table ..."; }
			$sth = $dbh->prepare("insert into MappingStats (libraryid, totalreads, mappedreads, alignmentrate, deletions, insertions, junctions, date ) values (?,?,?,?,?,?,?,?)");
			$sth ->execute($libraryidnumber, $totalreads, $mapped, $alignrate, $deletions, $insertions, $junctions, $date) or die "\nERROR:\t Complication in MappingStats table, consult documentation\n";
			if ($mappingtool) { print " Done\n"; }
			#metadata table
			if ($mappingtool) { print "NOTICE:\t Importing $mappingtool alignment information for $libraryidnumber to Metadata table ..."; }
			$sth = $dbh->prepare("insert into TheMetadata (libraryid, refgenome, annfile, stranded, sequences, mappingtool,status ) values (?,?,?,?,?,?,?)");
			$sth ->execute($libraryidnumber, $refgenomename, $annotationfile, $stranded,$sequences, $mappingtool,"done") or die "\nERROR:\t Complication in TheMetadata table, consult documentation\n";
			
			#Insert DataSyntaxes
			if ($mparameters) {
				$sth = $dbh->prepare("insert into Syntaxes (libraryid, mapsyntax ) values (?,?)");
				$mparameters =~ s/\"//g;
				$sth ->execute($libraryidnumber, $mparameters) or die "\nERROR:\t Complication in Syntaxes table, consult documentation\n";
			}
			if ($mappingtool) { print " Done\n"; }
			
			GENES_FPKM($libraryidnumber);
			READ_COUNT($libraryidnumber);
			DBVARIANT($variantfile, $libraryidnumber);
			print " Done\n";
					#variant annotation specifications
			if ($vep) {
				print "TASK:\t Importing Variant annotation from VEP => $folder2import\n"; #status
				print "NOTICE:\t Importing $libraryidnumber - Variant Annotation to VariantAnnotation table ...";
				VEPVARIANT($vepfile, $libraryidnumber); print " Done\n";
				NOSQL($libraryidnumber);
			}
			if ($annovar) {
				print "TASK:\t Importing Variant annotation from ANNOVAR => $folder2import\n"; #status
				print "NOTICE:\t Importing $libraryidnumber - Variant Annotation to VariantAnnotation table ...";
				ANNOVARIANT($annofile, $libraryidnumber); print " Done\n";
				NOSQL($libraryidnumber);
			}
		} else { #end unless found in MappingStats table
			print "NOTICE:\t $libraryidnumber already in MappingStats table... Moving on \n";
			$sth = $dbh->prepare("select libraryid from TheMetadata where libraryid = $libraryidnumber"); $sth->execute(); $found = $sth->fetch();
			unless ($found) {
				print "NOTICE:\t Importing $mappingtool alignment information for $libraryidnumber to Metadata table ...";
				$sth = $dbh->prepare("insert into TheMetadata (libraryid, refgenome, annfile, stranded, sequences, mappingtool,status) values (?,?,?,?,?,?,?)");
				$sth ->execute($libraryidnumber, $refgenomename, $annotationfile, $stranded,$sequences, $mappingtool,"done") or die "\nERROR:\t Complication in TheMetadata table, consult documentation\n";
				
				#Insert DataSyntaxes
				$sth = $dbh->prepare("insert into Syntaxes (libraryid, mapsyntax ) values (?,?)");
				$mparameters =~ s/\"//g;
				$sth ->execute($libraryidnumber, $mparameters) or die "\nERROR:\t Complication in Syntaxes table, consult documentation\n";
				print " Done\n";
				
			} #end else found in MappingStats table
			#toggle options
			GENES_FPKM($libraryidnumber); #GENES
			READ_COUNT($libraryidnumber); #READCOUNTS
			my $variantstatus = $dbh->selectrow_array("select status from VariantSummary where libraryid = $libraryidnumber and status = 'done'");
			unless ($variantstatus){ #checking if completed in VariantSummary table
				$verbose and print "NOTICE:\t Removed incomplete records for $libraryidnumber in all Variants tables\n";
				$sth = $dbh->prepare("delete from VariantAnnotation where libraryid = $libraryidnumber"); $sth->execute();
				$sth = $dbh->prepare("delete from VariantResult where libraryid = $libraryidnumber"); $sth->execute();
				$sth = $dbh->prepare("delete from VariantSummary where libraryid = $libraryidnumber"); $sth->execute();
				DBVARIANT($variantfile, $libraryidnumber);
				print " Done\n";
				#variant annotation specifications
				if ($vep) {
					print "TASK:\t Importing Variant annotation from VEP => $folder2import\n"; #status
					print "NOTICE:\t Importing $libraryidnumber - Variant Annotation to VariantAnnotation table ...";
					VEPVARIANT($vepfile, $libraryidnumber); print " Done\n";
					NOSQL($libraryidnumber);
				}
				if ($annovar) {
					print "TASK:\t Importing Variant annotation from ANNOVAR => $folder2import\n"; #status
					print "NOTICE:\t Importing $libraryidnumber - Variant Annotation to VariantAnnotation table ...";
					ANNOVARIANT($annofile, $libraryidnumber); print " Done\n";
					NOSQL($libraryidnumber);
				}
			} else {
				print "NOTICE:\t $libraryidnumber already in VariantResult table... Moving on \n";
				if ($vep || $annovar) {
					my $variantstatus = $dbh->selectrow_array("select annversion from VariantSummary where libraryid = $libraryidnumber");
					my $nosqlstatus = $dbh->selectrow_array("select nosql from VariantSummary where libraryid = $libraryidnumber");
					unless ($variantstatus && $nosqlstatus){ #if annversion or nosqlstatus is not specified
						$verbose and print "NOTICE:\t Removed incomplete records for $libraryidnumber in VariantAnnotation table\n";
						$sth = $dbh->prepare("delete from VariantAnnotation where libraryid = $libraryidnumber"); $sth->execute();
						if ($vep) {
							print "TASK:\t Importing Variant annotation from VEP => $folder2import\n"; #status
							print "NOTICE:\t Importing $libraryidnumber - Variant Annotation to VariantAnnotation table ...";
							VEPVARIANT($vepfile, $libraryidnumber); print " Done\n";
							NOSQL($libraryidnumber);
						}
						if ($annovar) {
							print "TASK:\t Importing Variant annotation from ANNOVAR => $folder2import\n"; #status
							print "NOTICE:\t Importing $libraryidnumber - Variant Annotation to VariantAnnotation table ...";
							ANNOVARIANT($annofile, $libraryidnumber); print " Done\n";
							NOSQL($libraryidnumber);
						}
					} else { #end unless annversion is previously specified
						print "NOTICE:\t $libraryidnumber already in VariantAnnotation table... Moving on\n";
					}
				} #end if annversion is previously specified
			} #end unless it's already in variants table
		} #end unless default option is specified
	} #unless & else exists in Mappingstats						
	else {
		pod2usage("FAILED: \"$libraryidnumber\" sample information is not in the database. Make sure the sample has been inputted from the website.");
	} #end if data in sample table							
								
	print ("SUCCESS: Import of RNA Seq analysis information in \"$folder2import\"\n");						
} 

if ($delete){ #delete section
	($ardea, $ibis) = fastbit();
	my (%KEYDELETE, $decision);
	my ($i,$alldelete) = (0,0);
	unless ($log) { print "JOB:\t Deleting Existing Records in Database\n"; } #status
	$dbh = mysql(); #connect to mysql
  	$sth = $dbh->prepare("select libraryid from BirdLibraries where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
	if ($found) {
		unless ($log) { print "NOTICE:\t This module deletes records from ALL database systems. Proceed with caution\n"; }
		$sth = $dbh->prepare("select libraryid from MappingStats where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
    		if ($found) {
			$i++; $KEYDELETE{$i} = "Alignment Information";
		}
		$sth = $dbh->prepare("select libraryid from GenesSummary where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
    		if ($found) {
			$i++; $KEYDELETE{$i} = "Expression Information";
		}
		$sth = $dbh->prepare("select libraryid from VariantSummary where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
    		if ($found) {
			$i++; $KEYDELETE{$i} = "Variant Information";
		}
		unless ($log) {
			print "--------------------------------------------------------------------------\n";
    		print "The following details match the libraryid '$delete' provided\n";
    		foreach (sort {$a <=> $b} keys %KEYDELETE) { print "  ", uc($_),"\.  $KEYDELETE{$_}\n";}
		}
		$KEYDELETE{0} = "ALL information relating to '$delete'";
		unless ($log) {
			print "  0\.  ALL information relating to '$delete'\n";
			print "--------------------------------------------------------------------------\n";
			print "Choose which information you want remove (multiple options separated by comma) or press ENTER to leave ? ";
			chomp ($decision = (<>)); print "\n";
		} else {$decision = $log; }
		if (length $decision >0) {
			my @allverdict = split(",",$decision);
			foreach my $verdict (sort {$b<=>$a} @allverdict) {
				if (exists $KEYDELETE{$verdict}) {
					print "NOTICE:\t Deleting $KEYDELETE{$verdict}\n";
					if ($verdict == 0) {$alldelete = 1;}
					if ($KEYDELETE{$verdict} =~ /^Variant/ || $alldelete == 1) {
						if ($alldelete == 1){
							if ($KEYDELETE{$i} =~ /^Variant/) { $i--;
								my $ffastbit = fastbit_path(); #connect to fastbit
								my $vfastbit = $ffastbit."/variant-information"; # specifying the variant section.
								print "NOTICE:\t Deleting records for $delete in Variant tables ";
								$sth = $dbh->prepare("delete from VariantAnnotation where libraryid = '$delete'"); $sth->execute(); print ".";
								$sth = $dbh->prepare("delete from VariantResult where libraryid = '$delete'"); $sth->execute(); print ".";
								$sth = $dbh->prepare("delete from VariantSummary where libraryid = '$delete'"); $sth->execute(); print ".";
								my $execute = "$ibis -d $vfastbit -y \"libraryid = '$delete'\" -z";
								`$execute`; print ".";
								`rm -rf $vfastbit/*sp $vfastbit/*old $vfastbit/*idx $vfastbit/*dic $vfastbit/*int `; #removing old indexes
								`$ibis -d $vfastbit -query "select genename, geneid, genetype, transcript, feature, codonchange, aminoacidchange, libraryid, chrom, tissue, species, consequence, existingvariant, source"`; #create a new index based on genename
								print " Done\n";
							}
						} else {
							my $ffastbit = fastbit_path(); #connect to fastbit
							my $vfastbit = $ffastbit."/variant-information"; # specifying the variant section.
							print "NOTICE:\t Deleting records for $delete in Variant tables ";
							$sth = $dbh->prepare("delete from VariantAnnotation where libraryid = '$delete'"); $sth->execute(); print ".";
							$sth = $dbh->prepare("delete from VariantResult where libraryid = '$delete'"); $sth->execute(); print ".";
							$sth = $dbh->prepare("delete from VariantSummary where libraryid = '$delete'"); $sth->execute(); print ".";
							my $execute = "$ibis -d $vfastbit -y \"libraryid = '$delete'\" -z";
							`$execute`; print ".";
							`rm -rf $vfastbit/*sp $vfastbit/*old $vfastbit/*idx $vfastbit/*dic $vfastbit/*int `; #removing old indexes
							`$ibis -d $vfastbit -query "select genename, geneid, genetype, transcript, feature, codonchange, aminoacidchange, libraryid, chrom, tissue, species, consequence, existingvariant, source"`; #create a new index
							print " Done\n";
						}
					}
					if ($KEYDELETE{$verdict} =~ /^Expression/ || $alldelete ==1 ) {
						if ($alldelete == 1){
							if ($KEYDELETE{$i} =~ /^Expression/) { $i--;
								my $ffastbit = fastbit_path();  #connect to fastbit
								my $gfastbit = $ffastbit."/gene-information"; # specifying the gene section.
								my $cfastbit = $ffastbit."/gene_count-information"; #specifying the gene count section
								
								print "NOTICE:\t Deleting records for $delete in Gene tables ";
								$sth = $dbh->prepare("delete from GenesSummary where libraryid = '$delete'"); $sth->execute(); print ".";
							
								#deleting gene-information from fastbit
								my $execute = "$ibis -d -v $gfastbit -y \"libraryid = '$delete'\" -z";
								`$execute`; print ".";
								`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
								`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
								
								#deleting gene_counts information from fastbit
								$execute = "$ibis -d -v $cfastbit -y \"libraryid = '$delete'\" -z";
								`$execute`; print ".";
								`rm -rf $cfastbit/*sp $cfastbit/*old $cfastbit/*idx $cfastbit/*dic $cfastbit/*int `; #removing old indexes
								`$ibis -d $cfastbit -query "select genename, libraryid, tissue, species"`; #create a new index based on genename
								
								print " Done\n";
							}
						} else {
							my $ffastbit = fastbit_path();  #connect to fastbit
							my $gfastbit = $ffastbit."/gene-information"; # specifying the gene section.
							my $cfastbit = $ffastbit."/gene_count-information"; #specifying the gene count section
								
							print "NOTICE:\t Deleting records for $delete in Gene tables ";
							$sth = $dbh->prepare("delete from GenesSummary where libraryid = '$delete'"); $sth->execute(); print ".";
							
							#deleting gene-information from fastbit
							my $execute = "$ibis -d -v $gfastbit -y \"libraryid = '$delete'\" -z";
							`$execute`; print ".";
							`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
							`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
								
							#deleting gene_counts information from fastbit
							$execute = "$ibis -d -v $cfastbit -y \"libraryid = '$delete'\" -z";
							`$execute`; print ".";
							`rm -rf $cfastbit/*sp $cfastbit/*old $cfastbit/*idx $cfastbit/*dic $cfastbit/*int `; #removing old indexes
							`$ibis -d $cfastbit -query "select genename, libraryid, tissue, species"`; #create a new index based on genename
								
							print " Done\n";
						}
					}
					if ($KEYDELETE{$verdict} =~ /^Alignment/ || $alldelete ==1 ) {
						if ($alldelete == 1){
							if ($KEYDELETE{$i} =~ /^Alignment/) { $i--;
								$sth = $dbh->prepare("select libraryid from GenesSummary where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
								unless ($found) {
									$sth = $dbh->prepare("select libraryid from VariantSummary where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
									unless ($found) {
										print "NOTICE:\t Deleting records for $delete in Mapping tables .";
										$sth = $dbh->prepare("delete from TheMetadata where libraryid = '$delete'"); $sth->execute(); print ".";
										$sth = $dbh->prepare("delete from Syntaxes where libraryid = '$delete'"); $sth->execute(); print ".";
										$sth = $dbh->prepare("delete from MappingStats where libraryid = '$delete'"); $sth->execute();  print ".";
										print " Done\n";
									} else { print "ERROR:\t Variant Information relating to '$delete' is in the database. Delete Variant Information first\n";}
								} else { print "ERROR:\t Expression Information Relating to '$delete' still present in the database. Delete Expression Information first\n";}
							}
						} else {
							$sth = $dbh->prepare("select libraryid from GenesSummary where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
							unless ($found) {
								$sth = $dbh->prepare("select libraryid from VariantSummary where libraryid = '$delete'"); $sth->execute(); $found = $sth->fetch();
								unless ($found) {
									print "NOTICE:\t Deleting records for $delete in Mapping tables .";
									$sth = $dbh->prepare("delete from TheMetadata where libraryid = '$delete'"); $sth->execute(); print ".";
									$sth = $dbh->prepare("delete from Syntaxes where libraryid = '$delete'"); $sth->execute(); print ".";
									$sth = $dbh->prepare("delete from MappingStats where libraryid = '$delete'"); $sth->execute();  print ".";
									print " Done\n";
								} else { print "ERROR:\t Variant Information relating to '$delete' is in the database. Delete Variant Information first\n";}
							} else { print "ERROR:\t Expression Information Relating to '$delete' still present in the database. Delete Expression Information first\n";}
						}
					}
				} else { print "ERROR:\t $verdict is an INVALID OPTION\n"; }
			}
		} else {print "NOTICE:\t No Option selected\n";}
	} else {
		print "NOTICE:\t Information relating to '$delete' is not in the database. Good-bye\n";
	}
}# -----------------------------------
# CREATING EMAIL NOTIFICATION
if ($sendemail) { NOTIFICATION("Completed Job",0,0); }
# -----------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - -S U B R O U T I N E S- - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub processArguments {
	my @commandline = @ARGV;
	GetOptions ( "h|help" => \$help, "man|manual" => \$manual, "del|delete=s" => \$delete, "send|email" => \$sendemail) or pod2usage ();

	$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
  	$manual and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);  
	unless ($delete) {
		@ARGV==1 or pod2usage("ERROR:\t No file specified\n"); 
		if ($ARGV[0]) { $folder2import = $ARGV[0]; } #library to transfer
	} # end unless delete file is specified

	#setup log file
	$nosql = @{ open_unique(".nosqlimport.txt") }[1]; `rm -rf $nosql`;
	$vnosql = @{ open_unique(".nosqlvimport.txt") }[1]; `rm -rf $vnosql`;
	$gnosql = @{ open_unique(".nosqlgimport.txt") }[1]; `rm -rf $gnosql`;
	$cnosql = @{ open_unique(".nosqlcimport.txt") }[1]; `rm -rf $cnosql`;
}


sub LOGFILE { #subroutine for getting metadata
	if ($fastqcfolder) { #if the fastqcfolder exist
		my ($fastqcfilename, $parentfastqc);
		if ($fastqcfolder =~ /zip$/) { #making sure if it's a zipped file
			`unzip $fastqcfolder`;
			$parentfastqc = fileparse($fastqcfolder, qr/\.[^.]*(\.zip)?$/);
			$fastqcfilename = fileparse($fastqcfolder, qr/\.[^.]*(\.zip)?$/)."/fastqc_data.txt";
		} else { $fastqcfilename = $fastqcfolder; } #else it will be the actual file
		$totalreads = `grep "Total Sequences" $fastqcfilename | awk -F" " '{print \$3}'`;
		if ($fastqcfolder =~ /zip$/) { `rm -rf $parentfastqc`; } #removing the unzipped folder
	}
	if ($bamfile){
		my $headerdetails = `samtools view -H $bamfile | grep -m 1 "\@PG" | head -1`;
		if ($headerdetails =~ /\sCL\:/) { #making sure mapping tool has the tool information or not 
			$headerdetails =~ /\@PG\s*ID\:(\S*).*VN\:(\S*)\s*CL\:(.*)/; $mappingtool = $1." v".$2; $mparameters = $3;
		} 
		else { $headerdetails =~ /\@PG\s*ID\:(\S*).*VN\:(\S*)/; $mappingtool = $1." v".$2; }
		if ($mappingtool =~ /hisat/i) {
			$mparameters =~ /\-x\s(\S+)\s/;
			$refgenome = $1; #reference genome name
			$refgenomename = (split('\/', $refgenome))[-1];
			if ($mparameters =~ /-1/){ #paired-end reads
				$mparameters =~ /\-1\s(\S+)\s-2\s(\S+)"$/;
				my @nseq = split(",",$1); my @pseq = split(",",$2);
				foreach (@nseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
				foreach (@pseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
				chop $sequences;
			}
			elsif ($mparameters =~ /-U/){ #single-end reads
				$mparameters =~ /\-U\s(\S+)"$/;
				my @nseq = split(",",$1);
				foreach (@nseq){ $sequences .= ( (split('\/', $_))[-1] ).",";}
				chop $sequences;
			} #end if toggle for sequences
			$stranded = undef;
			$annotationfile = undef;
		} # end if working with hisat.
		elsif ($mappingtool =~ /tophat/i) {
			my $annotation; undef %ALL; my $no = 0;
			my @newgeninfo = split('\s', $mparameters);
			my $number = 1;
			while ($number <= $#newgeninfo) {
				unless ($newgeninfo[$number] =~ /-no-coverage-search/){
					if ($newgeninfo[$number] =~ /^\-/){
						my $old = $number++;
						$ALL{$newgeninfo[$old]} = $newgeninfo[$number];
					} else {
						unless (exists $ALL{$no}){
							$ALL{$no} = $newgeninfo[$number];
							$no++;
						}
					}
				}
				$number++;
			}
			unless ((exists $ALL{"-G"}) || (exists $ALL{"--GTF"})) {
				$annotationfile = undef;
			} else {
				if (exists $ALL{"-G"}){ $annotation = $ALL{"-G"} ; } else { $annotation = $ALL{"--GTF"};}
				$annotationfile = uc ( (split('\.',((split("\/", $annotation))[-1])))[-1] ); #(annotation file)
			}
			unless (exists $ALL{"--library-type"}) { $stranded = undef; } else { $stranded = $ALL{"--library-type"}; }
		
			$refgenome = $ALL{0}; my $seq = $ALL{1}; my $otherseq = $ALL{2};
			$refgenomename = (split('\/', $ALL{0}))[-1]; 
			unless(length($otherseq)<1){ #sequences
				$sequences = ( ( split('\/', $seq) ) [-1]).",". ( ( split('\/', $otherseq) ) [-1]);
			} else {
				$sequences = ( ( split('\/', $seq) ) [-1]);
			} #end if seq
		} # end if working with tophat
		elsif ($mappingtool =~ /star/i) {
			my ($annotation, $otherseq); undef %ALL; my $no = 0; $mparameters =~ s/\s+/ /g;
			my @newgeninfo = split('\s', $mparameters);
			my $number = 1;
			while ($number <= $#newgeninfo) {
				unless ($newgeninfo[$number] =~ /-readFilesIn/){
					if ($newgeninfo[$number] =~ /^\-/){
						my $old = $number++;
						$ALL{$newgeninfo[$old]} = $newgeninfo[$number];
					}
				} else {
					#my $old = $number++;
					my $seq = $newgeninfo[++$number]; 
					my $new = $number+1;
					unless ($newgeninfo[$new] =~ /^\-\-/) {
						$otherseq = $newgeninfo[$new]; #if paired reads
						$sequences = ( ( split('\/', $seq) ) [-1]).",". ( ( split('\/', $otherseq) ) [-1]);
						$number++;
					} else {
						$sequences = ( ( split('\/', $seq) ) [-1]);
					}
				} #working with sequence names
				$number++;
			}
			$annotationfile = undef;
			$stranded = undef;
			$refgenome = $ALL{"--genomeDir"};
			$refgenomename = (split('\/', $ALL{"--genomeDir"}))[-1];
		} # end if working with star
		else {
			$annotationfile = undef;
			$stranded = undef; $sequences = undef;
		}
	} # end if bamfile
	if ($kallistologfile){
		my $versionnumber = `cat $kallistologfile | grep "kallisto_version" | awk -F'":' '{print \$2}'`; $versionnumber =~ s/\"//g; $versionnumber = substr($versionnumber,0,-2);
		$diffexpress = "kallisto v$versionnumber";
		$gparameters = `cat $kallistologfile | grep "call" | awk -F'":' '{print \$2}'`; $gparameters =~ s/\"//g;
		$gparameters =~ /-i\s(\S+)/; $refgenome = $1; $refgenomename = (split('\/', $refgenome))[-1]; $annotationfile = $refgenomename;
		my @newgeninfo = split(/\s/, $gparameters);
		if ($gparameters =~ /--single/) {	
			$sequences = ( ( split('\/', ($newgeninfo[-1])) ) [-1]);
		} else { #paired end reads
			$sequences = ( ( split('\/', $newgeninfo[-1]) ) [-2]).",". ( ( split('\/', $newgeninfo[-1]) ) [-1]);
		}
		$stranded = undef; $sequences = undef; $mappingtool = undef;
	} # end if kallistologfile
	if ($salmonlogfile){
		my $versionnumber = `cat $salmonlogfile | grep "salmon_version" | awk -F'":' '{print \$2}'`; $versionnumber =~ s/\"//g; $versionnumber = substr($versionnumber,0,-2);
		$diffexpress = "salmon v$versionnumber";
		$gparameters = `cat $salmonlogfile`;
		$refgenome = `cat $salmonlogfile | grep "index" | awk -F'":' '{print \$2}'`; $refgenome =~ s/\"//g; $refgenome = substr($refgenome,0,-2);
		$refgenomename = (split('\/', $refgenome))[-1]; $annotationfile = $refgenomename;
		$sequences = `cat $salmonlogfile | grep "mate" | awk -F'":' '{print \$2}'`; $sequences =~ s/\"//g; $sequences = substr($sequences,0,-2);
		my @newgeninfo = split(/\n/, $sequences); $sequences = undef;
		foreach (@newgeninfo) { $sequences .= (( split('\/', ($_)) ) [-1]); }
		$stranded = undef; $sequences = undef; $mappingtool = undef;
	} #end if salmonlogfile
}

sub READ_COUNT { #subroutine for read counts
	#INSERT INTO DATABASE: #ReadCounts table
	($ardea, $ibis) = fastbit();
	if ($readcountfile || $starcountfile) {
		my $countcolumn = "-1";
		my ($countpreamble, $checkforpreamble, $readcount) = (0,0,0);
		$sth = $dbh->prepare("select countstatus from GenesSummary where libraryid = $_[0] and countstatus ='done'"); $sth->execute(); $found = $sth->fetch();
		unless ($found) {
			my $ffastbit = fastbit_path();
			my $cfastbit = $ffastbit."/gene_count-information"; # specifying the gene section.
			`$ibis -d $cfastbit -v -q 'select count(libraryid) where libraryid = $_[0]' -o $nosql`;
			open(IN,"<",$nosql);
			no warnings;
			chomp($readcount = <IN>);
			close (IN); `rm -rf $nosql`;
			
			#get the species name and the tissue from the database
			my $species = $dbh->selectrow_array("select species from BirdLibraries where libraryid = $_[0]");
			my $tissue = $dbh->selectrow_array("select tissue from BirdLibraries where libraryid = $_[0]");			
				
			#get type of input
			if ($readcountfile) {
				open(READ, "<", $readcountfile) or die "\nERROR:\t Can not open file $readcountfile\n";
			} else {
				open(READ, "<", $starcountfile) or die "\nERROR:\t Can not open file $starcountfile\n";
				$countcolumn = "1"; $checkforpreamble = "1";
			}
			open (NOSQL, ">$cnosql");
			while (<READ>) {
				chomp;
				my @allidgene = split("\t");
				my ($idgene, $idcount) = ($allidgene[0], $allidgene[$countcolumn]); 
				if ($idgene =~ /^[a-zA-Z0-9][a-zA-Z0-9]/) {
					if ($countpreamble == 0 || $countpreamble == 2) {
						$checkforpreamble = 1;
					}
					if ($checkforpreamble == 1) {
						print NOSQL "'$idgene','$species','$tissue',$idcount,$_[0]\n"; #to fastbit
					} 
				} else {
					if ($_ =~ /featurecounts/i) {
						my ($program, $command) = split(';');
						$counttool = (split(':', $program))[1];
						$cparameters = (split(':', $command))[1];
						$cparameters =~ s/\"//g;
					} elsif ($_ =~ /^\_/) {
						$counttool = "htseqcount";
						$cparameters = undef;
					} elsif ($_ =~ /^N\_/) {
						$counttool = "STAR quantMode";
						$cparameters = undef;
					}
				}
				$countpreamble++;
			} close (READ);
			close NOSQL; #end of nosql portion
			
			if ($readcount < $countpreamble && $readcount != 0) {
				$verbose and print "NOTICE:\t Removed incomplete records for $_[0] in ReadCounts table\n";
				`$ibis -d $cfastbit -v -y \"libraryid = $_[0]\" -z`;
				`rm -rf $cfastbit/*sp $cfastbit/*old $cfastbit/*idx $cfastbit/*dic $cfastbit/*int `; #removing old indexes
			}	
			
			print "NOTICE:\t Importing $counttool raw counts information for $_[0] to ReadCounts table ...";
			#import into ReadCounts table;
			($ardea, $ibis) = fastbit();
			my $ffastbit = fastbit_path(); #connect to fastbit path
			my $cfastbit = $ffastbit."/gene_count-information"; # specifying the gene section.
			my $execute = "$ardea -d $cfastbit -m 'genename:text,species:text,tissue:text,readcount:double,libraryid:int' -t $cnosql";
			`$execute` or die "\nERROR\t: Complication importing RawCounts information to FastBit, contact $AUTHOR\n";
			`rm -rf $cfastbit/*sp`; #removeing old indexes
			`$ibis -d $cfastbit -query "select genename, libraryid, tissue, species"`; #create a new index based on genename
			`chmod 777 $cfastbit && rm -rf $cnosql`;
			
			$sth = $dbh->prepare("update GenesSummary set countstool = '$counttool', countstatus = 'done' where libraryid= $_[0]"); $sth ->execute(); #updating GenesSummary table.
			$sth = $dbh->prepare("update Syntaxes set countsyntax = '$cparameters' where libraryid= $_[0]"); $sth ->execute(); #updating Syntaxes table.
			print " Done \n";	
				
		} else { #found and completed
			print "NOTICE:\t $_[0] already in ReadCounts table... Moving on \n";
		}
	}
}		
		
sub GENES_FPKM { #subroutine for getting gene information
	($ardea, $ibis) = fastbit();
	#INSERT INTO DATABASE: #GenesSummary table
	my $ffastbit = fastbit_path();
	my $gfastbit = $ffastbit."/gene-information"; # specifying the gene section.
	
	$sth = $dbh->prepare("select libraryid from GenesSummary where libraryid = $_[0]"); $sth->execute(); $found = $sth->fetch();
	unless ($found) { 
		print "NOTICE:\t Importing $_[0] to GenesSummary table\n";
		$sth = $dbh->prepare("insert into GenesSummary (libraryid,date) values (?,?)");
		$sth ->execute($_[0], $date) or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
	} else {
		print "NOTICE:\t $_[0] already in GenesSummary table... Moving on \n";
	}
	my $genecount = 0;
	$sth = $dbh->prepare("select genestatus from GenesSummary where libraryid = $_[0] and genestatus ='done'"); $sth->execute(); $found = $sth->fetch();
	unless ($found) {
		`$ibis -d $gfastbit -v -q 'select count(libraryid) where libraryid = $_[0]' -o $nosql`;
		open(IN,"<",$nosql);
		no warnings;
		chomp($genecount = <IN>);
		close (IN); `rm -rf $nosql`;
		
		#get the species name and the tissue from the database
		my $species = $dbh->selectrow_array("select species from BirdLibraries where libraryid = $_[0]");
		my $tissue = $dbh->selectrow_array("select tissue from BirdLibraries where libraryid = $_[0]");
		
		#working on the individual fles
		if ($genesfile){ #working with genes.fpkm_tracking file
			#cufflinks expression tool name
			$diffexpress = "Cufflinks";
			$genes = `cat $genesfile | wc -l`; if ($genes >=2){ $genes--;} else {$genes = 0;} #count the number of genes
			$sth = $dbh->prepare("update GenesSummary set genes = $genes, diffexpresstool = '$diffexpress' where libraryid= $_[0]"); $sth ->execute(); #updating GenesSummary table.
			unless ($genes == $genecount) {
				unless ($genecount == 0 ) {
					$verbose and print "NOTICE:\t Removed incomplete records for $_[0] in Genes-NoSQL\n";
		      `$ibis -d $gfastbit -v -y \"libraryid = $delete\" -z`;
					`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
				}
				print "NOTICE:\t Importing $diffexpress expression information for $_[0] to Genes-NoSQL ...";
				#import into FPKM table;
				open(FPKM, "<", $genesfile) or die "\nERROR:\t Can not open file $genesfile\n";
				open (NOSQL, ">$gnosql");
				while (<FPKM>){
					chomp;
					my ($track, $class, $ref_id, $gene, $gene_name, $tss, $locus, $length, $coverage, $fpkm, $fpkm_low, $fpkm_high, $fpkm_stat ) = split /\t/;
					unless ($track eq "tracking_id"){ #check & specifying undefined variables to null
						if($coverage =~ /-/){$coverage = 0;} if (length $gene < 1) { $gene = "NULL"; } if (length $gene_name < 1) {$gene_name = "NULL";}
						my ($chrom_no, $chrom_start, $chrom_stop) = $locus =~ /^(.+)\:(.+)\-(.+)$/; $chrom_start++;

						print NOSQL "'$_[0]','$chrom_no','$gene','$gene_name','$species','$fpkm_stat','$tissue',$coverage,0,$fpkm,$fpkm_low,$fpkm_high,$chrom_start,$chrom_stop\n";
					}
				} close FPKM;
				close NOSQL; #end of nosql portion
		
				my $execute = "$ardea -d $gfastbit -m 'chrom:text,geneid:text,genename:text,species:text,fpkmstatus:char,tissue:text,coverage:double,tpm:double,fpkm:double,fpkmconflow:double,fpkmconfhigh:double,libraryid:int,chromstart:int,chromstop:int' -t $gnosql";
				`$execute` or die "\nERROR\t: Complication importing Expression information to FastBit, contact $AUTHOR\n";
				`rm -rf $gfastbit/*sp`; #removeing old indexes
				`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
				`chmod 777 $gfastbit && rm -rf $gnosql`;
				
				print " Done\n";
				#set GenesSummary to Done
				$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
				$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
			} else {
					$verbose and print "NOTICE:\t $_[0] already in Genes-NoSQL ... Moving on \n";
					$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
					$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
			}
		} elsif ($transcriptsgtf){ #using gtf files
			#differential expression tool names
			if (`head -n 1 $transcriptsgtf` =~ /cufflinks/i) { #working with cufflinks transcripts.gtf file
				$diffexpress = "Cufflinks";
				open(FPKM, "<", $transcriptsgtf) or die "\nERROR:\t Can not open file $transcriptsgtf\n";
				(%ARFPKM,%CHFPKM, %BEFPKM, %CFPKM, %DFPKM, %DHFPKM, %DLFPKM, %cfpkm, %dfpkm, %dhfpkm, %dlfpkm)= ();
				my $i=1;
				while (<FPKM>){
					chomp;
					my ($chrom_no, $tool, $typeid, $chrom_start, $chrom_stop, $qual, $orn, $idk, $therest ) = split /\t/;
					if ($typeid =~ /^transcript/){ #check to make sure only transcripts are inputed
						my %Drest = ();
						foreach (split("\";", $therest)) { $_ =~ s/\s+|\s+//g;my($a, $b) = split /\"/; $Drest{$a} = $b;}
						my $dstax;
						if (length $Drest{'gene_id'} > 1) {
							$dstax = "$Drest{'gene_id'}-$chrom_no";} else {$dstax = "xxx".$i++."-$chrom_no";}
						if (exists $CHFPKM{$dstax}){ #chromsome stop
							if ($chrom_stop > $CHFPKM{$dstax}) {
								$CHFPKM{$dstax} = $chrom_stop;
							}
						}else {
							$CHFPKM{$dstax} = $chrom_stop;
						}
						if (exists $BEFPKM{$dstax}){ #chromsome start
							if ($chrom_start < $BEFPKM{$dstax}) {
								$BEFPKM{$dstax} = $chrom_start;
							}
						}else {
							$BEFPKM{$dstax} = $chrom_start;
						}
						unless (exists $CFPKM{$dstax}{$Drest{'cov'}}){ #coverage
							$CFPKM{$dstax}{$Drest{'cov'}}= $Drest{'cov'};
						}unless (exists $DFPKM{$dstax}{$Drest{'FPKM'}}){ #FPKM
							$DFPKM{$dstax}{$Drest{'FPKM'}}= $Drest{'FPKM'};
						}
						unless (exists $DHFPKM{$dstax}{$Drest{'conf_hi'}}){ #FPKM_hi
							$DHFPKM{$dstax}{$Drest{'conf_hi'}}= $Drest{'conf_hi'};
						}
						unless (exists $DLFPKM{$dstax}{$Drest{'conf_lo'}}){ #FPKM_lo
							$DLFPKM{$dstax}{$Drest{'conf_lo'}}= $Drest{'conf_lo'};
						}
						$ARFPKM{$dstax}= "'$_[0]','$chrom_no','$Drest{'gene_id'}','$Drest{'gene_id'}'";
					}
				} close FPKM;
				#sorting the fpkm values and coverage results.
				foreach my $a (keys %DFPKM){
					my $total = 0;
					foreach my $b (keys %{$DFPKM{$a}}) { $total = $b+$total; }
					$dfpkm{$a} = $total;
				}
				foreach my $a (keys %CFPKM){
					my $total = 0;
					foreach my $b (keys %{$CFPKM{$a}}) { $total = $b+$total; }
					$cfpkm{$a} = $total;
				}
				foreach my $a (keys %DHFPKM){
					my $total = 0;
					foreach my $b (keys %{$DHFPKM{$a}}) { $total = $b+$total; }
					$dhfpkm{$a} = $total;
				}
				foreach my $a (keys %DLFPKM){
					my $total = 0;
					foreach my $b (keys %{$DLFPKM{$a}}) { $total = $b+$total; }
					$dlfpkm{$a} = $total;
				}
				#end of sort.
				#insert into database.
				$genes = scalar (keys %ARFPKM);
				$sth = $dbh->prepare("update GenesSummary set genes = $genes, diffexpresstool = '$diffexpress' where libraryid= $_[0]"); $sth ->execute(); #updating GenesSummary table.
				unless ($genes == $genecount) {
					unless ($genecount == 0 ) {
						$verbose and print "NOTICE:\t Removed incomplete records for $_[0] in Genes-NoSQL \n";
						`$ibis -d $gfastbit -v -y \"libraryid = $_[0]\" -z`;
						`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
					}
					print "NOTICE:\t Importing $diffexpress expression information for $_[0] to Genes-NoSQL ...";
					#import into FPKM table;
					open (NOSQL, ">$gnosql");
					foreach my $a (keys %ARFPKM){
						print NOSQL "$ARFPKM{$a},'$species','NULL','$tissue',$cfpkm{$a},0,$dfpkm{$a},$dlfpkm{$a},$dhfpkm{$a},$BEFPKM{$a},$CHFPKM{$a}\n";
					}
					close NOSQL; #end of nosql portion
					
					my $execute = "$ardea -d $gfastbit -m 'chrom:text,geneid:text,genename:text,species:text,fpkmstatus:char,tissue:text,coverage:double,tpm:double,fpkm:double,fpkmconflow:double,fpkmconfhigh:double,libraryid:int,chromstart:int,chromstop:int' -t $gnosql";
					`$execute` or die "\nERROR\t: Complication importing Expression information to FastBit, contact $AUTHOR\n";
					`rm -rf $gfastbit/*sp`; #removeing old indexes
					`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
					`chmod 777 $gfastbit && rm -rf $gnosql`;
				
					print " Done\n";
					#set GenesSummary to Done
					$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
					$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	else {
						$verbose and print "NOTICE:\t $_[0] already in Genes-NoSQL... Moving on \n";	
						$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
						$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	
			}
			elsif (`head -n 1 $transcriptsgtf` =~ /stringtie/i) { #working with stringtie output
				$gparameters = substr( `head -n 1 $transcriptsgtf`,2,-1 );
				$diffexpress = substr( `head -n 2 $transcriptsgtf | tail -1`,2,-1 );
				open(FPKM, "<", $transcriptsgtf) or die "\nERROR:\t Can not open file $transcriptsgtf\n";
				(%ARFPKM,%CHFPKM, %BEFPKM, %CFPKM, %DFPKM, %TPM, %cfpkm, %dfpkm, %tpm)= ();
				my $i=1;
				while (<FPKM>){
					chomp;
					my ($chrom_no, $tool, $typeid, $chrom_start, $chrom_stop, $qual, $orn, $idk, $therest ) = split /\t/;
					if ($typeid && $typeid =~ /^transcript/){ #check to make sure only transcripts are inputed
						my %Drest = ();
						foreach (split("\";", $therest)) { $_ =~ s/\s+|\s+//g;my($a, $b) = split /\"/; $Drest{$a} = $b;}
						my $dstax;
						if (length $Drest{'gene_id'} > 1) {
							$dstax = "$Drest{'gene_id'}-$chrom_no";} else {$dstax = "xxx".$i++."-$chrom_no";}
						if (exists $CHFPKM{$dstax}){ #chromsome stop
							if ($chrom_stop > $CHFPKM{$dstax}) {
								$CHFPKM{$dstax} = $chrom_stop;
							}
						}else {
							$CHFPKM{$dstax} = $chrom_stop;
						}
						if (exists $BEFPKM{$dstax}){ #chromsome start
							if ($chrom_start < $BEFPKM{$dstax}) {
								$BEFPKM{$dstax} = $chrom_start;
							}
						}else {
							$BEFPKM{$dstax} = $chrom_start;
						}
						unless (exists $CFPKM{$dstax}{$Drest{'cov'}}){ #coverage
							$CFPKM{$dstax}{$Drest{'cov'}}= $Drest{'cov'};
						}unless (exists $DFPKM{$dstax}{$Drest{'FPKM'}}){ #FPKM
							$DFPKM{$dstax}{$Drest{'FPKM'}}= $Drest{'FPKM'};
						}
						unless (exists $TPM{$dstax}{$Drest{'TPM'}}){ #TPM
							$TPM{$dstax}{$Drest{'TPM'}}= $Drest{'TPM'};
						}
						unless ($Drest{'ref_gene_name'}){
							$ARFPKM{$dstax}= "'$chrom_no','$Drest{'gene_id'}','NULL'";
						} else {
							$ARFPKM{$dstax}= "'$chrom_no','$Drest{'gene_id'}','$Drest{'ref_gene_name'}'";
						}
					}
				} close FPKM;
				#sorting the fpkm values and coverage results.
				foreach my $a (keys %DFPKM){
					my $total = 0;
					foreach my $b (keys %{$DFPKM{$a}}) { $total = $b+$total; }
					$dfpkm{$a} = $total;
				}
				foreach my $a (keys %CFPKM){
					my $total = 0;
					foreach my $b (keys %{$CFPKM{$a}}) { $total = $b+$total; }
					$cfpkm{$a} = $total;
				}
				foreach my $a (keys %TPM){
					my $total = 0;
					foreach my $b (keys %{$TPM{$a}}) { $total = $b+$total; }
					$tpm{$a} = $total;
				}
				#end of sort.
				#insert into database.
				$genes = scalar (keys %ARFPKM);
				$sth = $dbh->prepare("update GenesSummary set genes = $genes, diffexpresstool = '$diffexpress' where libraryid= $_[0]"); $sth ->execute(); #updating GenesSummary table.
				$gparameters =~ s/\"//g;
				$sth = $dbh->prepare("update Syntaxes set expsyntax = '$gparameters' where libraryid= $_[0]"); $sth ->execute(); #updating Syntaxes table.
			
				unless ($genes == $genecount) {
					unless ($genecount == 0 ) {
						$verbose and print "NOTICE:\t Removed incomplete records for $_[0] in Genes-NoSQL\n";
						`$ibis -d $gfastbit -v -y \"libraryid = $_[0]\" -z`;
						`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
					}
					print "NOTICE:\t Importing StringTie expression information for $_[0] to Genes-NoSQL ...";
					#import into FPKM table;
					open (NOSQL, ">$gnosql");
					foreach my $a (keys %ARFPKM){
						print NOSQL "$ARFPKM{$a},'$species','NULL','$tissue',$cfpkm{$a},$tpm{$a},$dfpkm{$a},0,0,$_[0],$BEFPKM{$a},$CHFPKM{$a}\n";
					}
					close NOSQL; #end of nosql portion
					
					my $execute = "$ardea -d $gfastbit -m 'chrom:text,geneid:text,genename:text,species:text,fpkmstatus:char,tissue:text,coverage:double,tpm:double,fpkm:double,fpkmconflow:double,fpkmconfhigh:double,libraryid:int,chromstart:int,chromstop:int' -t $gnosql";
					`$execute` or die "\nERROR\t: Complication importing Expression information to FastBit, contact $AUTHOR\n";
					`rm -rf $gfastbit/*sp`; #removing old indexes
					`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
					`chmod 777 $gfastbit && rm -rf $gnosql`;
				
					print " Done\n";
					#set GenesSummary to Done
					$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
					$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	else {
						$verbose and print "NOTICE:\t $_[0] already in Genes-NoSQL ... Moving on \n";
						$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
						$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	
			}	else {
				die "\nFAILED:\tCan not identify source of Genes Expression File '$transcriptsgtf', consult documentation.\n";
			}
		} elsif ($diffexpress =~ /kallisto/i) { #working with kallisto output
				open(FPKM, "<", $kallistofile) or die "\nERROR:\t Can not open file $kallistofile\n";
				(%TPM, %tpm) = ();
				my $i=1;
				while (<FPKM>){
					chomp;
					my ($targetid, $length, $eff, $est, $tpm ) = split /\t/;
					$TPM{$targetid}{$tpm} = $tpm;
				} close FPKM;
				#sorting the fpkm values and coverage results.
				foreach my $a (keys %TPM){
					my $total = 0;
					foreach my $b (keys %{$TPM{$a}}) { $total = $b+$total; }
					$tpm{$a} = $total;
				}
				#end of sort.
				#insert into database.
				$genes = scalar (keys %TPM);
				$sth = $dbh->prepare("update GenesSummary set genes = $genes, diffexpresstool = '$diffexpress' where libraryid= $_[0]"); $sth ->execute(); #updating GenesSummary table.
				$gparameters =~ s/\"//g; 
				$sth = $dbh->prepare("insert into Syntaxes (libraryid, expressionsyntax ) values (?,?)");
				$sth ->execute($_[0], $gparameters) or die "\nERROR:\t Complication in Syntaxes table, consult documentation\n";
				unless ($genes == $genecount) {
					unless ($genecount == 0 ) {
						$verbose and print "NOTICE:\t Removed incomplete records for $_[0] in Genes-NoSQL\n";
						`$ibis -d $gfastbit -v -y \"libraryid = $_[0]\" -z`;
						`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
					}
					print "NOTICE:\t Importing Kallisto expression information for $_[0] to Genes-NoSQL ...";
					#import into FPKM table;
					open (NOSQL, ">$gnosql");
					foreach my $a (keys %TPM){
						print NOSQL "'NULL','$a','$a','$species','NULL','$tissue',0,$tpm{$a},0,0,0,$_[0],0,0\n";
					}
					close NOSQL; #end of nosql portion
					
					my $execute = "$ardea -d $gfastbit -m 'chrom:text,geneid:text,genename:text,species:text,fpkmstatus:char,tissue:text,coverage:double,tpm:double,fpkm:double,fpkmconflow:double,fpkmconfhigh:double,libraryid:int,chromstart:int,chromstop:int' -t $gnosql";
					`$execute` or die "\nERROR\t: Complication importing Expression information to FastBit, contact $AUTHOR\n";
					`rm -rf $gfastbit/*sp`; #removing old indexes
					`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
					`chmod 777 $gfastbit && rm -rf $gnosql`;
				
					print " Done\n";
					#set GenesSummary to Done
					$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
					$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	else {
						$verbose and print "NOTICE:\t $_[0] already in Genes-NoSQL ... Moving on \n";
						$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
						$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	
		}	elsif ($diffexpress =~ /salmon/i) { #working with salmon output
				open(FPKM, "<", $salmonfile) or die "\nERROR:\t Can not open file $salmonfile\n";
				(%TPM, %tpm) = ();
				my $i=1;
				while (<FPKM>){
					chomp;
					my ($targetid, $length, $eff, $tpm, $est ) = split /\t/;
					$TPM{$targetid}{$tpm} = $tpm;
				} close FPKM;
				#sorting the fpkm values and coverage results.
				foreach my $a (keys %TPM){
					my $total = 0;
					foreach my $b (keys %{$TPM{$a}}) { $total = $b+$total; }
					$tpm{$a} = $total;
				}
				#end of sort.
				#insert into database.
				$genes = scalar (keys %TPM);
				$sth = $dbh->prepare("update GenesSummary set genes = $genes, diffexpresstool = '$diffexpress' where libraryid= $_[0]"); $sth ->execute(); #updating GenesSummary table.
				$gparameters =~ s/\"//g;
				$sth = $dbh->prepare("insert into Syntaxes (libraryid, expressionsyntax ) values (?,?)");
				$sth ->execute($_[0], $gparameters) or die "\nERROR:\t Complication in Syntaxes table, consult documentation\n";
				
				unless ($genes == $genecount) {
					unless ($genecount == 0 ) {
						$verbose and print "NOTICE:\t Removed incomplete records for $_[0] in Genes-NoSQL \n";
						`$ibis -d $gfastbit -v -y \"libraryid = $_[0]\" -z`;
						`rm -rf $gfastbit/*sp $gfastbit/*old $gfastbit/*idx $gfastbit/*dic $gfastbit/*int `; #removing old indexes
					}
					print "NOTICE:\t Importing Salmon expression information for $_[0] to Genes-NoSQL ...";
					#import into FPKM table;
					open (NOSQL, ">$gnosql");
					foreach my $a (keys %TPM){
						print NOSQL "'$_[0]','NULL','$a','$a','$species','NULL','$tissue',0,$tpm{$a},0,0,0,0,0\n";
					}
					close NOSQL; #end of nosql portion
					
					my $execute = "$ardea -d $gfastbit -m 'chrom:text,geneid:text,genename:text,species:text,fpkmstatus:char,tissue:text,coverage:double,tpm:double,fpkm:double,fpkmconflow:double,fpkmconfhigh:double,libraryid:int,chromstart:int,chromstop:int' -t $gnosql";
					`$execute` or die "\nERROR\t: Complication importing Expression information to FastBit, contact $AUTHOR\n";
					`rm -rf $gfastbit/*sp`; #removing old indexes
					`$ibis -d $gfastbit -query "select genename, geneid, libraryid, chrom, tissue, species"`; #create a new index based on genename
					`chmod 777 $gfastbit && rm -rf $gnosql`;
				
					print " Done\n";
					#set GenesSummary to Done
					$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
					$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	else {
						$verbose and print "NOTICE:\t $_[0] already in Genes-NoSQL ... Moving on \n";
						$sth = $dbh->prepare("update GenesSummary set genestatus = 'done' where libraryid = $_[0]");
						$sth ->execute() or die "\nERROR:\t Complication in GenesSummary table, consult documentation\n";
				}	
		}
	} else {
		$verbose and print "NOTICE:\t $_[0] already in Genes-NoSQL ... Moving on \n";
		
	}
}

sub DBVARIANT {
	my $toolvariant; undef $vparameters;
	unless($_[0]) { print ("ERROR:\t Can not find variant file. make sure variant file with suffix '.vcf' is present\nNOTICE:\t Moving on from importing variants to TransAtlasDB ..."); }
	else { open(VARVCF,$_[0]) or die ("\nERROR:\t Can not open variant file $_[0]\n");  
		while (<VARVCF>) {
			chomp;
			if (/^\#/) {
				unless (/\#INFO/ || /\#FORMAT/ || /\#contig/ || /#FORMAT/) {
					if (/Version/) {
						if ($_ =~ /GATK/) {
							$_ =~ /ID\=(.*)\,.*Version\=(.*)\,Date.*CommandLineOptions="(.*)">$/;
							$toolvariant = "GATK v.$2,$1";
							$varianttool = "GATK";
							$vparameters = $3;
						} elsif ($_ =~ /samtools/) {
							$_ =~ /Version\=(.*)\+./;
							$toolvariant = "samtools v.$1";
							$varianttool = "samtools";
						}
					} elsif (/Command/) {
						$_ =~ /Command=(.*)$/;
						unless ($vparameters) {
							$vparameters = $1;
						} else {
							$vparameters .= " | $1";
						}
					} #end assigning toolvariant
				}
			} else {
				my @chrdetails = split "\t";
				my @morechrsplit = split(';', $chrdetails[7]);
				if (((split(':', $chrdetails[9]))[0]) eq '0/1'){$verd = "heterozygous";}
				elsif (((split(':', $chrdetails[9]))[0]) eq '1/1'){$verd = "homozygous";}
				elsif (((split(':', $chrdetails[9]))[0]) eq '1/2'){$verd = "heterozygous alternate";}
				$VCFhash{$chrdetails[0]}{$chrdetails[1]} = "$chrdetails[3]|$chrdetails[4]|$chrdetails[5]|$verd";
			}
		} close VARVCF;
		$sth = $dbh->prepare("insert into VariantSummary ( libraryid, varianttool, date) values (?,?,?)");
		$sth ->execute($_[1], $toolvariant, $date) or die "\nERROR:\t Complication in VariantSummary table, consult documentation\n";;
		$vparameters =~ s/\"//g;
		$sth = $dbh->prepare("update Syntaxes set variantsyntax = '$vparameters' where libraryid = $_[1]");
		$sth ->execute();
	
		#VARIANT_RESULTS
		print "NOTICE:\t Importing $varianttool variant information for $_[1] to VariantResult table ...";
		
		#converting to threads
		undef %HASHDBVARIANT; my $ii = 0;
		foreach my $abc (sort keys %VCFhash) {
			foreach my $def (sort {$a <=> $b} keys %{ $VCFhash{$abc} }) {
				my @vcf = split('\|', $VCFhash{$abc}{$def});
				if ($vcf[3] =~ /,/){
					my $first = split(",",$vcf[1]);
					if (length $vcf[0] == length $first){ $itvariants++; $itsnp++; $variantclass = "SNV"; }
					elsif (length $vcf[0] < length $first) { $itvariants++; $itindel++; $variantclass = "insertion"; }
					else { $itvariants++; $itindel++; $variantclass = "deletion"; }
				}
				elsif (length $vcf[0] == length $vcf[1]){ $itvariants++; $itsnp++; $variantclass = "SNV"; }
				elsif (length $vcf[0] < length $vcf[1]) { $itvariants++; $itindel++; $variantclass = "insertion"; }
				else { $itvariants++; $itindel++; $variantclass = "deletion"; }
			
				#putting variants info into a hash table
				my @hashdbvariant = ($_[1], $abc, $def, $vcf[0], $vcf[1], $vcf[2], $variantclass, $vcf[3]); 
				$HASHDBVARIANT{$ii++} = [@hashdbvariant];
			}
		}
			
		#update variantsummary with counts
		$sth = $dbh->prepare("update VariantSummary set totalvariants = $itvariants, totalsnps = $itsnp, totalindels = $itindel where libraryid= $_[1]");
		$sth ->execute();
	
		my @hashdetails = keys %HASHDBVARIANT; #print "First $#hashdetails\n"; die;
		undef @VAR; undef @threads;
		push @VAR, [ splice @hashdetails, 0, 200 ] while @hashdetails; #sub the files to multiple subs
		$queue = new Thread::Queue();
		my $builder=threads->create(\&main); #create thread for each subarray into a thread 
		push @threads, threads->create(\&dbvarprocessor) for 1..5; #execute 5 threads
		$builder->join; #join threads
		foreach (@threads){$_->join;}
	
		#update variantsummary with counts
		$sth = $dbh->prepare("update VariantSummary set status = 'done' where libraryid= $_[1]");
		$sth ->execute();
		$sth->finish();
	}
}
sub dbvarprocessor { my $query; while ($query = $queue->dequeue()){ dbvarinput(@$query); } }
sub dbvarinput {
	foreach my $a (@_) { 
		$dbh = mysql(); #connect to mysql
		$sth = $dbh->prepare("insert into VariantResult ( libraryid, chrom, position, refallele, altallele, quality, variantclass, zygosity ) values (?,?,?,?,?,?,?,?)");
		$sth -> execute(@{$HASHDBVARIANT{$a}}) or die "\nERROR:\t Complication in VariantResult table, consult documentation\n";
	}
}

sub VEPVARIANT {
	my ($chrom, $position);
	if($_[0]){ open(VEP,$_[0]) or die ("\nERROR:\t Can not open vep file $_[0]\n"); } else { die ("\nERROR:\t Can not find VEP file. make sure vep file with suffix '.vep.txt' is present\n"); }
	my $ii = 0; undef %HASHDBVARIANT;
	while (<VEP>) {
		chomp;
		unless (/^\#/) {
			unless (/within_non_coding_gene/i || /coding_unknown/i) {
				my @veparray = split "\t"; #14 columns
				my @extraarray = split(";", $veparray[13]);
				foreach (@extraarray) { my @earray = split "\="; $extra{$earray[0]}=$earray[1]; }
				unless (length($veparray[0]) >1) {
					my @indentation = split(":", $veparray[1]);
					$chrom = $indentation[0]; $position = $indentation[1];
					unless ( $extra{'VARIANT_CLASS'} =~ /SNV/i || $extra{'VARIANT_CLASS'} =~ /substitution/i ){
						if($position =~ /\-/) { $position = (split("\-", $position))[0]; }
						unless ($extra{'VARIANT_CLASS'} =~ /insertion/) { $position--; }
						if ($extra{'VARIANT_CLASS'} =~ /indel/){ #check if the chromosomal number matches the VCF file
							my $check = 0;
							$sth = $dbh->prepare("select chrom, position from VariantResult where libraryid = $_[1] and chrom = '$chrom' and position = $position"); $sth->execute(); $found = $sth->fetch();
							unless($found) { $position++; }
						}
					}
				} else {
					my @indentation = split("_", $veparray[0]);
					if ($#indentation > 2) { $chrom = $indentation[0]."_".$indentation[1]; $position = $indentation[2]; }
					else { $chrom = $indentation[0]; $position = $indentation[1]; }
					unless ( $extra{'VARIANT_CLASS'} =~ "SNV" or $extra{'VARIANT_CLASS'} =~ "substitution" ){ $position--; }
					else {
						my @poly = split("/",$indentation[$#indentation]);
						unless ($#poly > 1){ unless (length ($poly[0]) == length($poly[1])){ $position--; } }
					}
				}
				my $geneid = $veparray[3];
				my $transcriptid = $veparray[4];
				my $featuretype = $veparray[5];
				my $consequence = $veparray[6]; 
				if ($consequence =~ /NON_(.*)$/){ $consequence = "NON".$1; } elsif ($consequence =~ /STOP_(.*)$/) {$consequence = "STOP".$1; }
				my $pposition = $veparray[9];
				my $aminoacid = $veparray[10];
				my $codons = $veparray[11];
				my $dbsnp = $veparray[12];
				my $locate = "$_[1],$chrom,$position,$consequence,$geneid,$pposition";
				if ( exists $VEPhash{$locate} ) {
					unless ( $VEPhash{$locate} eq $locate ){ die "\nERROR:\t Duplicate annotation in VEP file, consult documentation\n"; }
				} else {
					$VEPhash{$locate} = $locate;
					if (exists $extra{'SYMBOL'}) { $extra{'SYMBOL'} = uc($extra{'SYMBOL'}); }
					my @hashdbvep = ($_[1], $chrom, $position, $extra{'VARIANT_CLASS'}, $consequence, $geneid, $pposition, $dbsnp, $extra{'SOURCE'}, $extra{'SYMBOL'}, $transcriptid, $featuretype, $extra{'BIOTYPE'}, $aminoacid, $codons);
					$HASHDBVARIANT{$ii++} = [@hashdbvep];
					$DBSNP{$chrom}{$position} = $dbsnp; #updating dbsnp	
				}
			}
		} else { if (/API (version \d+)/){ $annversion = $1;} } #getting VEP version
	}
	close VEP;
	
	#update convert to threads
	my @hashdetails = keys %HASHDBVARIANT; #print "First $#hashdetails\n"; die;
	undef @VAR; undef @threads;
	push @VAR, [ splice @hashdetails, 0, 200 ] while @hashdetails; #sub the files to multiple subs
	$queue = new Thread::Queue();
	my $builder=threads->create(\&main); #create thread for each subarray into a thread 
	push @threads, threads->create(\&dbveprocessor) for 1..5; #execute 5 threads
	$builder->join; #join threads
	foreach (@threads){$_->join;}
	`cat $vnosql* >.temp$vnosql; rm -rf $vnosql*; mv .temp$vnosql $vnosql;`; #merging all files
	
	$sth = $dbh->prepare("update VariantSummary set annversion = 'VEP $annversion' where libraryid = $_[1]"); $sth ->execute();
}
sub dbveprocessor { my $query; while ($query = $queue->dequeue()){ dbvepinput(@$query); } }
sub dbvepinput {
	$vcount++;
	my $newvcount = $vnosql."$vcount";
	foreach my $a (@_) { 
		$dbh = mysql(); #connect to mysql
		$sth = $dbh->prepare("insert into VariantAnnotation ( libraryid, chrom, position, consequence, geneid, proteinposition, source, genename, transcript, feature, genetype, aminoacidchange, codonchange ) values (?,?,?,?,?,?,?,?,?,?,?,?,?)");
		my ($vsample,$vchrom, $vposition, $vclass, $vconsequence, $vgeneid, $vpposition, $vdbsnp) = @{$HASHDBVARIANT{$a}}[0..7];
		my @vrest = @{$HASHDBVARIANT{$a}}[8..$#{$HASHDBVARIANT{$a}}];
		$sth->execute($vsample,$vchrom, $vposition, $vconsequence, $vgeneid, $vpposition, @vrest) or die "\nERROR:\t Complication in VariantAnnotation table, consult documentation\n";
		$sth = $dbh->prepare("update VariantResult set variantclass = '$vclass' where libraryid = $vsample and chrom = '$vchrom' and position = $vposition"); $sth ->execute() or die "\nERROR:\t Complication in updating VariantResult table, consult documentation\n";
		
		#DBSNP
		if (exists $DBSNP{$vchrom}{$vposition}) {
			$sth = $dbh->prepare("update VariantResult set existingvariant = '$vdbsnp' where libraryid = $vsample and chrom = '$vchrom' and position = $vposition"); $sth ->execute();
			delete $DBSNP{$vchrom}{$vposition};
		}
		#NOSQL portion
		@nosqlrow = $dbh->selectrow_array("select * from vw_variantnosql where libraryid = $vsample and chrom = '$vchrom' and position = $vposition and consequence = '$vconsequence' and geneid = '$vgeneid' and proteinposition = '$vpposition'");
		$showcase = undef; 
		foreach my $variables (0..$#nosqlrow){
			if ($variables == 2) { $nosqlrow[$variables] = $vdbsnp; }
			if (!($nosqlrow[$variables]) ||(length($nosqlrow[$variables]) < 1) || ($nosqlrow[$variables] =~ /^\-$/) ){
				$nosqlrow[$variables] = "NULL";
			}
			if ($variables < 17) {
				$nosqlrow[$variables] =~ s/^'|'$//g;
				$showcase .= "'$nosqlrow[$variables]',";
			}
			else {
				$showcase .= "$nosqlrow[$variables],";
			}
		}
		chop $showcase; $showcase .= "\n";
		open (NOSQL, ">>$newvcount"); print NOSQL $showcase; close NOSQL; #end of nosql portion
		undef %extra; #making sure extra is blank
	}
}

sub ANNOVARIANT {
	my (%REFGENE, %ENSGENE, %CONTENT);
	my $ii = 0; undef %HASHDBVARIANT;
	if($_[0]){ open(ANNOVAR,$_[0]) or die ("\nERROR:\t Can not open annovar file $_[0]\n"); } else { die ("\nERROR:\t Can not find annovar file. make sure annovar file with suffix '.multianno.txt' is present\n"); }
	my @annocontent = <ANNOVAR>; close ANNOVAR; 
	my @header = split("\t", lc($annocontent[0]));

	#getting headers
	foreach my $no (5..$#header-1){ #checking if ens and ref is specified
		my @tobeheader = split('\.', $header[$no]);
		if ($tobeheader[1] =~ /refgene/i){
			$REFGENE{$tobeheader[0]} = $header[$no];
		} elsif ($tobeheader[1] =~ /ensgene/i){ 
			$ENSGENE{$tobeheader[0]} = $header[$no];
		} else {
			die "ERROR:\t Do not understand notation '$tobeheader[1]' provided. Contact $AUTHOR \n";
		}
	} #end foreach dictionary
	#convert content to a hash array.
	my $counter = $#header+1;
	@annocontent = @annocontent[1..$#annocontent];
	foreach my $rowno (0..$#annocontent) {
		unless ($annocontent[$rowno] =~ /intergenic.*NONE,NONE/i) {
			my @arrayrow = split("\t", $annocontent[$rowno], $counter);
			foreach my $colno (0..$#header) {
				$CONTENT{$rowno}{$header[$colno]} = $arrayrow[$colno];
			}
			$CONTENT{$rowno}{'position'} = (split("\t",$arrayrow[$#arrayrow]))[4];
		} # end unless var annotation is useless
	} #end getting column position o entire file
	#working with ENS
	if (exists $ENSGENE{'func'}) {
		foreach my $newno (sort {$a<=>$b} keys %CONTENT){
			my $pposition = "-"; my $consequence = ""; my $transcript = ""; my $aminoacid = "-"; my $codons = "-";
			if ($CONTENT{$newno}{$ENSGENE{'func'}} =~ /^exonic/i) {
				$consequence = $CONTENT{$newno}{$ENSGENE{'exonicfunc'}};
				unless ($consequence =~ /unknown/i){         
					my @acontent = split(",", $CONTENT{$newno}{$ENSGENE{'aminoacidchange'}});
					my @ocontent = split (":", $acontent[$#acontent]);
					$transcript = $ocontent[1] =~ s/\.\d//g;
					foreach (@ocontent){
						if (/^c\.[a-zA-Z]/) {
							my ($a, $b, $c) = $_ =~ /^c\.(\S)(\d+)(\S)$/; 
							$codons = $a."/".$c;
						} elsif (/^c\.[0-9]/) {
							$codons = $_;
						} elsif (/^p\.\S\d+\S$/){ 
							my ($a, $b, $c) = $_ =~ /^p\.(\S)(\d+)(\S)$/;
							$pposition = $b;
							if ($a eq $c) {$aminoacid = $a;} else {$aminoacid = $a."/".$b;}
						} elsif (/^p\.\S\d+\S+$/) {
							my $a = $_ =~ /^p\.\S(\d+)\S+$/;
							$pposition = $a;
							$aminoacid = $_;
						}
					} #end foreach @ocontent
				} else {next;} 
			} else {
				$consequence = $CONTENT{$newno}{$ENSGENE{'func'}};
			}
			unless ($consequence =~ /ncRNA/i) {
				$consequence = uc($consequence);
				if ($consequence eq "UTR5"){ $consequence = "5PRIME_UTR";}
				if ($consequence eq "UTR3"){ $consequence = "3PRIME_UTR";} 
			}
			$CONTENT{$newno}{$ENSGENE{'gene'}} =~ s/\.\d//g; 
			my $locate = "$_[1],$CONTENT{$newno}{'chr'}, $CONTENT{$newno}{'position'},$consequence,$CONTENT{$newno}{$ENSGENE{'gene'}},$pposition";
			if ( exists $ANNOhash{$locate} ) {
				unless ( $ANNOhash{$locate} eq $locate ){ die "\nERROR:\t Duplicate annotation in ANNOVAR file, contact $AUTHOR\n"; }
			} else {
				$ANNOhash{$locate} = $locate;
				my @hashdbanno = ($_[1], $CONTENT{$newno}{'chr'}, $CONTENT{$newno}{'position'}, $consequence, $CONTENT{$newno}{$ENSGENE{'gene'}}, $pposition, 'Ensembl', $transcript, $aminoacid, $codons);
				$HASHDBVARIANT{$ii++} = [@hashdbanno];
			} #end if annohash locate
		} # end foreach looking at content
	} #end if ENSGENE

	#working with REF
	if (exists $REFGENE{'func'}) {
		foreach my $newno (sort {$a<=>$b} keys %CONTENT){
			my $pposition = "-"; my $consequence = ""; my $transcript = ""; my $aminoacid = ""; my $codons = "";
			if ($CONTENT{$newno}{$REFGENE{'func'}} =~ /^exonic/i) {
				$consequence = $CONTENT{$newno}{$REFGENE{'exonicfunc'}};
				unless ($consequence =~ /unknown/i){         
					my @acontent = split(",", $CONTENT{$newno}{$REFGENE{'aminoacidchange'}});
					my @ocontent = split (":", $acontent[$#acontent]);
					$transcript = $ocontent[1] =~ s/\.\d//g;
					foreach (@ocontent){
						if (/^c\.[a-zA-Z]/) {
							my ($a, $b, $c) = $_ =~ /^c\.(\S)(\d+)(\S)$/; 
							$codons = $a."/".$c;
						} elsif (/^c\.[0-9]/) {
							$codons = $_;
						} elsif (/^p\.\S\d+\S$/){ 
							my ($a, $b, $c) = $_ =~ /^p\.(\S)(\d+)(\S)$/;
							$pposition = $b;
							if ($a eq $c) {$aminoacid = $a;} else {$aminoacid = $a."/".$b;}
						} elsif (/^p\.\S\d+\S+$/) {
							my $a = $_ =~ /^p\.\S(\d+)\S+$/;
							$pposition = $a;
							$aminoacid = $_;
						}
					}
				} else { next; }
			} else {
				$consequence = $CONTENT{$newno}{$REFGENE{'func'}};
			}
			unless ($consequence =~ /ncRNA/i) {
				$consequence = uc($consequence);
				if ($consequence eq "UTR5"){ $consequence = "5PRIME_UTR";}
				if ($consequence eq "UTR3"){ $consequence = "3PRIME_UTR";} 
			}
			$CONTENT{$newno}{$REFGENE{'gene'}} =~ s/\.\d//g;
			my $locate = "$_[1],$CONTENT{$newno}{'chr'}, $CONTENT{$newno}{'position'},$consequence,$CONTENT{$newno}{$REFGENE{'gene'}},$pposition";
			if ( exists $ANNOhash{$locate} ) {
				unless ( $ANNOhash{$locate} eq $locate ){ die "\nERROR:\t Duplicate annotation in ANNOVAR file, contact $AUTHOR\n"; }
			} else {
				$ANNOhash{$locate} = $locate;
				if (exists $CONTENT{$newno}{$REFGENE{'gene'}}) { $CONTENT{$newno}{$REFGENE{'gene'}} = uc($CONTENT{$newno}{$REFGENE{'gene'}}); }
				my @hashdbanno = ($_[1], $CONTENT{$newno}{'chr'}, $CONTENT{$newno}{'position'}, $consequence, $CONTENT{$newno}{$REFGENE{'gene'}}, $pposition, 'RefSeq', $CONTENT{$newno}{$REFGENE{'gene'}}, $transcript, $aminoacid, $codons);
				$HASHDBVARIANT{$ii++} = [@hashdbanno];
			} #end if annohash locate
		} # end foreach looking at content
	} #end if REFGENE
	
	#update convert to threads
	my @hashdetails = keys %HASHDBVARIANT; #print "First $#hashdetails\n"; die;
	undef @VAR; undef @threads;
	push @VAR, [ splice @hashdetails, 0, 200 ] while @hashdetails; #sub the files to multiple subs
	$queue = new Thread::Queue();
	my $builder=threads->create(\&main); #create thread for each subarray into a thread 
	push @threads, threads->create(\&dbannorocessor) for 1..5; #execute 5 threads
	$builder->join; #join threads
	foreach (@threads){$_->join;}
	`cat $vnosql* >.temp$vnosql; rm -rf $vnosql*; mv .temp$vnosql $vnosql;`; #merging all files
	
	$sth = $dbh->prepare("update VariantSummary set annversion = 'ANNOVAR' where libraryid = $_[1]"); $sth ->execute(); #update database annversion :  ANNOVAR
}
sub dbannorocessor { my $query; while ($query = $queue->dequeue()){ dbannoinput(@$query); } }
sub dbannoinput {
	$vcount++;
	my $newvcount = $vnosql."$vcount";
	foreach my $a (@_) { 
		$dbh = mysql(); #connect to mysql
		if ($#{$HASHDBVARIANT{$a}} == 9) { $sth = $dbh->prepare("insert into VariantAnnotation ( libraryid, chrom, position, consequence, geneid, proteinposition, source, transcript, aminoacidchange, codonchange ) values (?,?,?,?,?,?,?,?,?,?)"); }
		else { $sth = $dbh->prepare("insert into VariantAnnotation ( libraryid, chrom, position, consequence, geneid, proteinposition, source, genename, transcript, aminoacidchange, codonchange ) values (?,?,?,?,?,?,?,?,?,?,?)"); }
		#print $#{$HASHDBVARIANT{$a}},"\n";
		my ($vsample,$vchrom, $vposition, $vconsequence, $vgeneid, $vpposition) = @{$HASHDBVARIANT{$a}}[0..5];
		my @vrest = @{$HASHDBVARIANT{$a}}[6..$#{$HASHDBVARIANT{$a}}];
		$sth->execute($vsample,$vchrom, $vposition, $vconsequence, $vgeneid, $vpposition, @vrest) or die "\nERROR:\t Complication in VariantAnnotation table, consult documentation\n";

		#NOSQL portion
		@nosqlrow = $dbh->selectrow_array("select * from vw_variantnosql where libraryid = $vsample and chrom = '$vchrom' and position = $vposition and consequence = '$vconsequence' and geneid = '$vgeneid' and proteinposition = '$vpposition'");
		$showcase = undef;
		foreach my $variables (0..$#nosqlrow){
			if (!($nosqlrow[$variables]) ||(length($nosqlrow[$variables]) < 1) || ($nosqlrow[$variables] =~ /^\-$/) ){
				$nosqlrow[$variables] = "NULL";
			}
			if ($variables < 17) {
				$showcase .= "'$nosqlrow[$variables]',";
			}
			else {
				$showcase .= "$nosqlrow[$variables],";
			}
		}
		chop $showcase; $showcase .= "\n";
		open (NOSQL, ">>$newvcount"); print NOSQL $showcase; close NOSQL; #end of nosql portion
		undef %extra; #making sure extra is blank
	}
}

sub NOSQL {
	print "TASK:\t Importing Variant annotation for $_[0] to NoSQL platform\n"; #status
	my $ffastbit = fastbit_path();  #connect to fastbit
	my $vfastbit = $ffastbit."/variant-information"; # specifying the variant section.
	print "NOTICE:\t Importing $_[0] - Variant Annotation to NoSQL '$ffastbit' ...";
	my $execute = "$ardea -d $vfastbit -m 'variantclass:char,zygosity:char,existingvariant:text,source:text,consequence:text,geneid:text,genename:text,transcript:text,feature:text,genetype:text,refallele:char,altallele:char,tissue:text,chrom:text,aminoacidchange:text,codonchange:text,species:text,quality:double,libraryid:int,position:int,proteinposition:int' -t $vnosql";
	`$execute` or die "\nERROR\t: Complication importing to FastBit, contact $AUTHOR\n";
	`rm -rf $vfastbit/*sp`; #removing old indexes
	`$ibis -d $vfastbit -query "select genename, geneid, genetype, transcript, feature, codonchange, aminoacidchange, libraryid, chrom, tissue, species, consequence, existingvariant, source"`; #create a new index
	`chmod 777 $vfastbit && rm -rf $vnosql`;
	$sth = $dbh->prepare("update VariantSummary set nosql = 'done' where libraryid = $_[0]"); $sth ->execute(); #update database nosql : DONE
	
	#removing records from MySQL
	$sth = $dbh->prepare("delete from VariantAnnotation where libraryid = $_[0]"); $sth->execute();

	#declare done
	print " Done\n";
}

sub NOTIFICATION {
   my $notification = '/home/modupe/.LOG/note.txt';
   if ($_[1] == 1) {
		open (NOTE, ">$notification");
		print NOTE "Subject: Update notes : $jobid\n\nCompleted library\t$_[2]\n";
		system "sendmail $EMAIL < $notification"; close NOTE;
      system "rm -rf $notification";
   } #end if
   else {
      open (NOTE, ">$notification");
      print NOTE "Subject: ". $_[0] .": $jobid\n\nName of log files\n\t$std_out\n\t$std_err\n";
      system "sendmail $EMAIL < $notification"; close NOTE;
      system "rm -rf $notification";
  }
}

sub main {
  foreach my $count (0..$#VAR) {
		while(1) {
			if ($queue->pending() < 100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
	}
	foreach(1..5) { $queue-> enqueue(undef); }
}
#--------------------------------------------------------------------------------

=pod

=head1 NAME

$0 -- Comprehensive pipeline : Inputs frnakenstein results in the database : transcriptatlas

=head1 SYNOPSIS

 importtodb.pl [arguments] [--help] [--delete] [--manual] <directory of files>

 Optional arguments:
        -del, --delete                  delete library information
				-h, --help                      print help message
        -m, --man                       print complete documentation


=head1 DESCRIPTION

Accepts all folders from frnakenstein output.
 
=head1 OPTIONS

=over 3

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2018 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage
=cut



