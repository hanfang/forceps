#include <iostream>
#include <string>
#include <sstream>
#include "align.hh"
#include "util.hh"

#include <boost/algorithm/string.hpp> 
#include <list>

using namespace boost; 
using namespace std;

void readseq(const string & filename, string & seq)
{
 	FILE * f = xfopen(filename, "r");
 	string hdr;

	if (!Fasta_Read(f, seq, hdr))
  	{
 		cerr << "ERROR: Can't read fasta record from: " << filename << endl;
	}

	cerr << "Loaded " << hdr << " from " << filename << endl;
}


int main(int argc, char ** argv)
{
 	string USAGE = "Usage: marlin [options] ref.fa qry.fq\n";

	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);

	stringstream helptext;
	helptext <<
	USAGE
	<<
	"\n"
	"options:\n"
	"   -v    : be verbose\n"
	"   -e    : Don't penalize if the ends of the sequences don't align\n"
	"   -u    : Enable Wobble base pairing, i.e. G-U base pairing (if this flag is not off, a missing value will be reported)\n"
	"   -d    : Distance to either end of the kmer (default: -d 3)\n"
	"   -r    : Do the reverse complement alignment as well\n"
	"\n";

	if (argc == 1)
	{
		cerr << helptext.str();
		exit(0);
	}

	bool errflg = false;
  	int ch;

	int ENDFREE = 0;
	int VERBOSE = 0;
	int REVERSE = 0;
	int ALLOWGU = 0;
  	int DISTOEND = 3;

	optarg = NULL;

	while (!errflg && ((ch = getopt (argc, argv, "verud:h")) != EOF))
  	{
    	switch (ch)
    	{
      		case 'v': VERBOSE = 1; break;
      		case 'e': ENDFREE = 1; break;
      		case 'r': REVERSE = 1; break; 
      		case 'u': ALLOWGU = 1; break; 
      		case 'd': DISTOEND = atoi(optarg); break;
      		case 'h': errflg = 1;  break;

      		case '?':
        	fprintf (stderr, "Unrecognized option -%c\n", optopt);

    		default:
    		errflg = true;
    	}

    	if (errflg)
    	{
    		cout << helptext.str();
    		exit (EXIT_FAILURE);
    	}
  	}

	if (argc < optind+2)
  	{
    	cerr << USAGE;
    	exit(1);
	}

  	// Parse options
  	FILE * ref_fa = xfopen(argv[optind], "r");
  	FILE * qry_fa = xfopen(argv[optind+1], "r");

  	// Load the refernece sequence
  	string ref_hdr;
  	string R;

  	cout << "#chr" << "\t" << "start" << "\t" << "end" << "\t" << "query" << "\t" << "strand" << "\t" << "edit_distance" << "\t" << "mis_match" << "\t" << "insertion" << "\t" << "deletion" << "\t" << "back_to_back" << "\t" << "near_edge" << "\t" << "Wooble" << endl;

  	while (Fasta_Read(ref_fa, R, ref_hdr))
	{
    	
    	// cout << "Loaded ref: " << ref_hdr << " [" << R.length() << "bp]" << endl << endl;
    	// cout << "==================================================================" << endl;

    	// Now load the query sequences and align

    	string qry_hdr;
    	string Q;
    	string Q_rc;
      
    	rewind(qry_fa);

    	while (Fasta_Read(qry_fa, Q, qry_hdr))
		{
			// cout << "Aligning qry: " << qry_hdr << " [" << Q.length() << "bp]" << endl;
	  		// cout << "==================================================================" << endl;
	  
	  		// Align in the forward orientations
	  		string R_aln;
	  		string Q_aln;
	  
	  	  	// Initialize integer variables
	  		int edit_distance = 0;
	  		int mis_match     = 0;
	  		int insertion     = 0;
	  		int deletion      = 0;
	  		int GtoU          = 0;
	  		int dis_to_end    = DISTOEND;
	  		// Initialize binary variables
	  		bool back_to_back = false;
	  		bool fwd_two 	  = false;
      		bool fwd_three 	  = false;
	  		bool edge         = false;
	  		bool hdl_end      = false;
	  		bool hdl_mid      = false;

	  		// cout << "R:  " << R << " [" << R.length() << "bp]" << endl;
	  		// cout << "Q:  " << Q << " [" << Q.length() << "bp]" << endl;
	  
	  		int dis_to_fiveprime = global_align_aff(R, Q, R_aln, Q_aln, ENDFREE, ALLOWGU, VERBOSE); 
	  		// cout << endl;
	  
	  		cout << "R': " << R_aln << " [" << R_aln.length() << "bp]" << endl;
	  		cout << "Q': " << Q_aln << " [" << Q_aln.length() << "bp]" << endl;
	  
	  		cout << "    ";
	  
	  		for (unsigned int i = 0; i < R_aln.length(); i++)
	    	{
	  			if (ALLOWGU == 0)
	  			{
	  				if (R_aln[i] != Q_aln[i]) 
	  				{
	  					edit_distance++;
		  				if ( i < Q_aln.length()-1 )
		  				{
		  					if ( R_aln[i+1] != Q_aln[i+1] )
		  					{
		  						fwd_two = true; 
		  						if ( R_aln[i+2] != Q_aln[i+2] )
		  						{
		  							fwd_three = true; 
		  						}
		  						else
		  						{ } // do nothing, pass
		  					}
		  				}
		  			}
		  		}			
		  		else
		  		{			
	  				if ( (R_aln[i] != Q_aln[i]) && ( !(R_aln[i] == 'T' && Q_aln[i] == 'C') ) )  
					{
						edit_distance++;
		  				if (i < Q_aln.length()-1 )
		 	   			{ 
		 	   				if ( R_aln[i+1] != Q_aln[i+1] && ( !(R_aln[i+1] == 'T' && Q_aln[i+1] == 'C') ) ) 
							{
				  				fwd_two = true; 
				  				if (R_aln[i+2] != Q_aln[i+2] && ( !(R_aln[i+2] == 'T' && Q_aln[i+2] == 'C') ) ) 
				    			{ 
				    				fwd_three = true; 
				    			}
			 					else 
			 	   				{ } // do nothing, pass 
							} 
		  	  			} 
		  			}
		  		}
		  		
		  		// Decide if the mutations occur near either edge
				if (i <= dis_to_end-1 || i >=Q_aln.length()-dis_to_end ) 
				{ 
					hdl_end = true; 
				}
				if (i > dis_to_end-1  && i < Q_aln.length()-dis_to_end ) 
				{ 
					hdl_mid = true; 
				}
	    		// cout << Q_aln.length() << endl ;	
				// cout << Q_aln.length()-dis_to_end+1 <<  endl ;

	    		// Mark the mismatches, indels, and gtou within the alignment
	    		if (R_aln[i] == Q_aln[i]) 
	    		{ 
	    			cout << " "; 
	    		}
	    		else if (R_aln[i] == '-') { cout << "^"; insertion++; }
	      		else if (Q_aln[i] == '-') { cout << "v"; deletion++; }
	      		else if ( (R_aln[i] == 'T') && (Q_aln[i] == 'C') && ALLOWGU == 1 ) { cout << "U"; GtoU++; }
	      		else                      { cout << "X"; mis_match++; }
	    		
	  			
	    	}
				
	    	// Parse the coordinates and retrieve the mapping locations
	  		vector<string> handler;
	  		vector<string> handler_level2;
			string chrom;	
			int start;

	  		boost::split(handler, ref_hdr, is_any_of(":"));
			for (int k = 0; k < handler.size(); ++k) 
			{
    			chrom = handler[0];

    			boost::split(handler_level2, handler[1], is_any_of("-"));
    			for (int m = 0; m < handler_level2.size(); ++m) 
    			{
    				start = atoi( handler_level2[0].c_str() );
    			}
    		}

    		int start_mapping = start + dis_to_fiveprime ; //
    		int end_mapping   = start_mapping + Q_aln.length() ; //

    		/*
    		// Convert the integers to string
    		string start_mapping_str;
    		ostringstream convert_start;
    		convert_start << start_mapping;
    		start_mapping_str = convert_start.str();

    		string end_mapping_str;
    		ostringstream convert_end;
    		convert_end << end_mapping;
    		end_mapping_str = convert_end.str();
				
			// Combine the string to report the alignment location
			std::string mapping_loc = chrom + ":" + start_mapping_str + "-" + end_mapping_str ; // 
			*/ 

			// Report if the mutations occur near the edge
			if (hdl_end == true && hdl_mid == false) 
			{ 
				edge = true;
			}

			// Report if there are back-to-back mutations in the alignment
	  		if ( (fwd_two == true && edit_distance == 2) || (fwd_three == true && edit_distance == 3) )
	  		{
	  			back_to_back = true;
	  		}

	  		cout << endl ;

	  		if (ALLOWGU == 0)
	  		{
	  			cout << "chr" << chrom << "\t" << start_mapping << "\t" << end_mapping << "\t" << qry_hdr << "\tFwd\t" << edit_distance << "\t" << mis_match << "\t" << insertion << "\t" << deletion << "\t" << back_to_back << "\t" << edge << "\t" << "." << endl;
			}
			else
			{
				cout << "chr" << chrom << "\t" << start_mapping << "\t" << end_mapping << "\t" << qry_hdr << "\tFwd\t" << edit_distance << "\t" << mis_match << "\t" << insertion << "\t" << deletion << "\t" << back_to_back << "\t" << edge << "\t" << GtoU << endl;
			}

		// Now align with reverse complement of Q	  
		if (REVERSE)
	  	{
	      	// Reset integer variables
	      	edit_distance = 0;
	      	mis_match     = 0;
	      	insertion     = 0;
	      	deletion      = 0;
	      	GtoU          = 0;
	      	// Reset binary variables
	      	back_to_back = false;
	      	fwd_two      = false;
	      	fwd_three    = false;
	      	edge         = false;
	      	hdl_end      = false;
	      	hdl_mid      = false;
	      	

	     	string Q_rc = rc_str(Q);

	      	// cout << "R:  " << R << " [" << R.length() << "bp]" << endl;
	      	// cout << "Q~: " << Q_rc << " [" << Q_rc.length() << "bp]" << endl;
	      
	      	int dis_to_fiveprime = global_align_aff(R, Q_rc, R_aln, Q_aln, ENDFREE, ALLOWGU, VERBOSE);
	     	// cout << endl;

	    	cout << "R':  " << R_aln << " [" << R_aln.length() << "bp]" << endl;
	      	cout << "Q~': " << Q_aln << " [" << Q_aln.length() << "bp]" << endl;
	      
	      	cout << "     ";
	      
	     	for (unsigned int i = 0; i < R_aln.length(); i++)
			{
		  		if (ALLOWGU == 0)
	  			{
	  				if (R_aln[i] != Q_aln[i]) 
	  				{
	  					edit_distance++;
		  				if ( i < Q_aln.length()-1 )
		  				{
		  					if ( R_aln[i+1] != Q_aln[i+1] )
		  					{
		  						fwd_two = true; 
		  						if ( R_aln[i+2] != Q_aln[i+2] )
		  						{
		  							fwd_three = true; 
		  						}
		  						else
		  						{ } // do nothing, pass
		  					}
		  				}
		  			}
		  		}			
		  		else
		  		{			
	  				if ( (R_aln[i] != Q_aln[i]) && ( !(R_aln[i] == 'T' && Q_aln[i] == 'C') ) )  
					{
						edit_distance++;
		  				if (i < Q_aln.length()-1 )
		 	   			{ 
		 	   				if ( R_aln[i+1] != Q_aln[i+1] && ( !(R_aln[i+1] == 'T' && Q_aln[i+1] == 'C') ) ) 
							{
				  				fwd_two = true; 
				  				if (R_aln[i+2] != Q_aln[i+2] && ( !(R_aln[i+2] == 'T' && Q_aln[i+2] == 'C') ) ) 
				    			{ 
				    				fwd_three = true; 
				    			}
			 					else 
			 	   				{ } // do nothing, pass 
							} 
		  	  			} 
		  			}
		  		}

		      	// Decide if the mutations occur near either edge
				if (i <= dis_to_end-1 || i >=Q_aln.length()-dis_to_end ) { hdl_end = true; }
				if (i > dis_to_end-1  && i < Q_aln.length()-dis_to_end ) { hdl_mid = true; }
		  		
		  
		  		if (R_aln[i] == Q_aln[i]) 
		  		{ 
		  			cout << " ";
		  		}
		  		else if (R_aln[i] == '-') { cout << "^"; insertion++;}
		  		else if (Q_aln[i] == '-') { cout << "v"; deletion++;}
		  		else if ( (R_aln[i] == 'T') && (Q_aln[i] == 'C') && ALLOWGU == 1 ) { cout << "U"; GtoU++; }
		  		else                      { cout << "X"; mis_match++;}
			}
	     
			// Parse the coordinates and retrieve the mapping locations
	  		vector<string> handler; //
	  		vector<string> handler_level2; //
			string chrom; //	  			
			int start; //

	  		boost::split(handler, ref_hdr, is_any_of(":"));  // 
			for (int k = 0; k < handler.size(); ++k) 
			{
    			chrom = handler[0];

    			boost::split(handler_level2, handler[1], is_any_of("-"));
    			for (int m = 0; m < handler_level2.size(); ++m) 
    			{
    				start = atoi( handler_level2[0].c_str() );
    				}
    		}

    		int start_mapping = start + dis_to_fiveprime ;
    		int end_mapping   = start_mapping + Q_aln.length() ;
    		/*
    		// Convert the integers to string
    		string start_mapping_str;
    		ostringstream convert_start;
    		convert_start << start_mapping;
    		start_mapping_str = convert_start.str();

    		string end_mapping_str;
    		ostringstream convert_end;
    		convert_end << end_mapping;
    		end_mapping_str = convert_end.str();
				
			// Combine the string to report the alignment location
			std::string mapping_loc = chrom + ":" + start_mapping_str + "-" + end_mapping_str ; 
			*/
			// Report if the mutations occur near the edge
	    	if (hdl_end == true && hdl_mid == false) 
	    	{
	    		edge = true;
	    	}
	      	if ( (fwd_two == true && edit_distance == 2) || (fwd_three == true && edit_distance == 3) )
	      	{
	      	back_to_back = true;
	      	}

	      	cout << endl ;	 
	      
			if (ALLOWGU == 0)
	  		{
	  			cout << "chr" << chrom << "\t" << start_mapping << "\t" << end_mapping << "\t" << qry_hdr << "\tRev\t" << edit_distance << "\t" << mis_match << "\t" << insertion << "\t" << deletion << "\t" << back_to_back << "\t" << edge << "\t" << "." << endl;
			}
			else
			{
				cout << "chr" << chrom << "\t" << start_mapping << "\t" << end_mapping << "\t" << qry_hdr << "\tRev\t" << edit_distance << "\t" << mis_match << "\t" << insertion << "\t" << deletion << "\t" << back_to_back << "\t" << edge << "\t" << GtoU << endl;
			}
	    }
	      	
	// cout << endl;
       
	}
	}

  return 0;	  
}
