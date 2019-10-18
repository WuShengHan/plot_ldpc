#!/usr/bin/perl -w

if( $#ARGV >= 0 ){
    print $ARGV[0];
    $infile = $ARGV[0];
}else{
    $infile = 'ecc_tb.v';
}
$output_file = "gen_$infile.pl";

open( in_f  , "$infile" )   || die " Not input file exist." ;
open( out_f , ">$output_file" ) || die " Not create output file." ;

print out_f "#!/usr/bin/perl -w\n\n";
printf out_f "open( OUT , \">tmp_%s\"); \n\n", $infile;
while(<in_f>){
    
    chomp;
    @char_str = split ( // , $_);
    print out_f "printf OUT \"";
    
    foreach $c_char (@char_str){
        if( $c_char eq '"' ){
            print out_f '\"';
        } elsif( $c_char eq '$' ){
            print out_f '\$';
        } elsif( $c_char eq '\\' ){
            print out_f '\\\\';
        } elsif( $c_char eq '@' ){
            print out_f '\@';
        } elsif( $c_char eq '%' ){
            print out_f '%%';
        } else {
            print out_f "$c_char";
        }    
    }
    print out_f "\\n\"; \n";
#    print out_f "printf OUT \" $_\\n \"; \n";
}
print out_f "\nclose(OUT);\n";

close(in_f)  || die " Can't close input file." ;
close(out_f) || die " Can't close output file.";