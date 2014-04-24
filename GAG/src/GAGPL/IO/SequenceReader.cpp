/********************************************************************
	created:	2012/06/24
	created:	24:6:2012   21:18
	filename: 	SequenceReader.cpp
	file path:	GAGPL\IO
	file base:	SequenceReader
	file ext:	cpp
	author:		Han Hu
	
	purpose:	
*********************************************************************/
#include <fstream>
#include <memory>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>

#include <GAGPL/IO/SequenceReader.h>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GLYCAN/GlycanComposition.h>
#include <GAGPL/MISC/Param.h>

namespace gag 
{
	namespace io
	{
		using namespace std;

		void SequenceReader::readSequenceStructure(const std::string& filename, GlycanSequencePtr seq)
		{
			std::ifstream seqfile;
			std::string ROW, flag;
			seqfile.open(filename.c_str());

			using namespace boost::xpressive;

			//std::auto_ptr<GlycanSequence> gs_pt;
			std::auto_ptr<Branch> bc_pt;
			
			GlycanComposition glycan_compo;
			std::string init_unit = "GlcA";

			sregex res_ori = '>' >> (s1= +_w);
			smatch what;

			//GlycanSequence seq;
			std::string construct_type = "SEQUENCE";

			if(seqfile.is_open()) {
				while(!seqfile.eof()) {

					getline(seqfile, ROW); // Read file by line.

					if(regex_search(ROW, what, res_ori)) {

						// flag stores the content of the headline for later use.
						flag = what[1];

						if(what[1] == "SEQUENCE") {  // multiple sequences can be appended in this format.
							
						} else if(what[1] == "END"){ // "END" represents the end of a glycan sequence.

						} else if(what[1] == "COMPOSITION") {
							// Construct the sequence by glycan code.
							construct_type = "CODE";
							//glycan_compo_pt = new GlycanComposition();

						} else if(what[1] == "MODIFICATION") {
							
						} else if(what[1] == "COMBINATION") {
							seq->buildByGAGComposition(glycan_compo, init_unit);
							construct_type = "CLEAR";
						} else if(what[1] == "INITIALIZATION") { 
							
						}else {
							std::cout << "Undefined headline!" << std::endl;
						}

					} else {

						if(flag == "SEQUENCE") { // Add the branches into sequence.	
							// The sequence is set as a backbone.
							// Split the sequence into modules.
							//sregex res = +_s;
							//sregex_token_iterator cur(ROW.begin(), ROW.end(), res, -1), end;
							
							//size_t mono_id = 0;
							//// Retrieve the branch id.
							//bc_pt = new Branch(boost::lexical_cast<size_t>(*cur));

							//for(++cur; cur != end; cur++)
							//{	
							//	std::string str = *cur;
							//	res = bos >> +_w >> eos;
							//	if(regex_match(str, what, res)) {
							//		Monosaccharide ms = mut.getMonosaccharideBySymbol(str);
							//		bc_pt->addUnit(ms);				
							//	} else {
							//		// If linkage.
							//		res = bos >> (s1= +_w) >> (s2= +_d) >> '-' >> (s3= +_d) >> eos;
							//		if(regex_match(str, what, res)){
							//			Linkage lk(mono_id, boost::lexical_cast<size_t>(what[2]), boost::lexical_cast<size_t>(what[3]), what[1]);
							//			bc_pt->addLinkage(lk);
							//			mono_id++;
							//		}
							//	}
							//}

						} else if(flag == "COMPOSITION") {
							sregex res = bos >> (s1= +_w) >> '\t' >> (s2= +_d);
							if(regex_match(ROW, what, res)){
								const std::string mono = what[1];
								int num = boost::lexical_cast<int>(what[2]);
								glycan_compo.addGlycanComposition(mono, num);
							}
							
						} else if(flag == "INITIALIZATION") {
							sregex res = bos >> (s1=+_w);
							if(regex_match(ROW, what, res)) init_unit = what[1];
						}  else if(flag == "MODIFICATION"){
							
								// s1 -- modification symbol
								// s2 -- branch id.
								// s3 -- monosaccharide id.
								// s4 -- site id.
								sregex res = bos >> (s1= +_w) >> '\t' >> (s2 = +_d) >> '\t' >>  (s3 = +_d) >> '\t' >> (s4 = +_d);
								if(regex_match(ROW, what, res)){
									std::string mod = what[1];
									ModificationPosition mod_pos(boost::lexical_cast<size_t>(what[2]), boost::lexical_cast<size_t>(what[3]), boost::lexical_cast<size_t>(what[4]));
									seq->addModification(what[1], mod_pos);
								}
						} else if(flag == "COMBINATION"){
							// s1 -- modification symbol
							// s2 -- modification number
							sregex res = bos >> (s1= +_w) >> '\t' >> (s2= +_d);
							if(regex_match(ROW, what, res)){
								std::string mod = what[1];
								int num = boost::lexical_cast<int>(what[2]);
								// Add modification number.
								seq->addModificationConstraint(mod, num);
							}

						}
					}

				} 
				
				seqfile.close();

			} else {
				throw std::runtime_error("Unable to open file!");
			}
			
			//return seq;
		}

    GlycanSequencePtr SequenceReader::readSequenceStructure(const std::string& filename, const std::string path)
    {
      GlycanSequencePtr gs = boost::make_shared<GlycanSequence>();
      return gs;
    }

		void SequenceReader::readSpectrum( const std::string& filename, RichList& spectrum )
		{
			std::ifstream infile(filename.c_str());

			std::string line;

			//RichList spec;
			if(infile.is_open())
			{
				while(std::getline(infile, line))
				{
					std::istringstream is;
					is.str(line);
					double k1, k2, k3, k4; 
					// The four columns are: m/z, intensity, resolution, s/n
					is >> k1 >> k2 >> k3 >> k4;
					//gag::Peak pk(k1, k2);
					RichPeakPtr pk = createRichPeak(k1, k2, k3, k4);
					spectrum.addPeak(pk);
				}
			}
			infile.close();
			infile.clear();	
		}

		void SequenceReader::readMonoPeakList( const std::string& filename, std::set<MonoPeakPtr>& mono_list )
		{
			std::ifstream infile(filename.c_str());

			std::string line;

			if(infile.is_open())
			{
				// Deal with title.
				std::getline(infile, line);	

				while(std::getline(infile, line))
				{
					std::istringstream is;
					is.str(line);
					// mz, intensity and charge state.
					double k1; double k2; int k3; 
					is >> k1 >> k2 >> k3;
					gag::MonoPeakPtr pk = boost::make_shared<MonoPeak>(k1, k2, k3);
					//double mass = msmath::calculateMass(k1, -1 * k2);
					//std::cout.precision(5);
					//std::cout << std::fixed << mass << std::endl;
					mono_list.insert(pk);
				}
			}
			infile.close();
			infile.clear();
		}

	//	void SequenceReader::writePredictionResults(const std::string& filename, std::string title, std::string seq_code, std::string m_code, const std::vector<double>& seq_scores)
	//	{
	//		static std::ofstream outfile(filename.c_str());
 //     static bool head = true;
 //     if(head) {
 //       if(outfile.is_open()) {
 //         outfile << "SEQ\tM_SEQ\tCoverage\tGP\tCost\tM_Coverage\tM_GP\tM_Cost\n";

 //       } else {
 //         throw std::runtime_error("Unable to open file!");
 //       }
 //       head = false;
 //     }

 //     outfile << seq_code << "\t" << m_code;

 //     for(auto iter = seq_scores.begin(); iter != seq_scores.end(); iter++)
 //     {
 //       outfile << 
 //     }

	//		if(outfile.is_open()) {
	//			// Title.
	//			outfile << "SEQ\tCoverage\tGP\tCost\tM_Coverage\tM_GP\tM_Cost\n";
	//			
	//			SequenceScoreList::const_iterator iter = seq_scores.begin();
	//			for(; iter != seq_scores.end(); iter++)
	//			{
	//				outfile << iter->first << "\t";
	//				const std::vector<double>& score_vec = iter->second;
	//				for(std::vector<double>::const_iterator score_it = score_vec.begin(); score_it != score_vec.end(); score_it++)
	//				{
	//					outfile << "\t" << *score_it;
	//				}
	//				outfile << "\n";
	//			}

	//			outfile.close();

	//		} 
	//	}

	//}
  }
	
}