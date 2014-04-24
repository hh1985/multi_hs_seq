#include <GAGPL/GAGLIBRARY/LibraryForest.h>
//#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
#include <GAGPL/FRAGMENTATION/Fragmentation.h>
#include <GAGPL/CHEMISTRY/Modifier.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <utility>
#include <boost/xpressive/xpressive.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
//#include <GAGPL/SPECTRUM/PseudoclusterFinderHelper.h>

using namespace std;
using namespace gag;

typedef std::multimap<double, int> MonoList;

void readSequence(std::vector<GlycanSequence>& gag_list, std::map<std::string, size_t>& mod_list)
{
	//MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	FunctionalGroupTable& fun_table = FunctionalGroupTable::Instance();
	ifstream seqfile;
	string ROW, flag;
	seqfile.open("../../test/Boons_11_structure.txt");

	using namespace boost::xpressive;
	GlycanSequence* gs_pt(new GlycanSequence);
	//typedef boost::shared_ptr<GlycanSequence> GlycanSequencePtr;
	//GlycanSequencePtr gs_pt = boost::make_shared<GlycanSequence>();
	Branch* bc_pt(new Branch);
	/*typedef boost::shared_ptr<Branch> BranchPtr;
	BranchPtr bc_pt = boost::make_shared<Branch>();*/

	sregex res_ori = '>' >> (s1= +_w);

	smatch what;

	if(seqfile.is_open()) {
		while(!seqfile.eof()) {
			getline(seqfile, ROW);

			if(regex_search(ROW, what, res_ori)) {
				std::string next_flag = what[1];
				//cout << flag << endl;
				if(next_flag == "END"){
					gs_pt->update();
					gag_list.push_back(*gs_pt);
					delete gs_pt;
				} else if(next_flag == "SEQUENCE" && !gag_list.empty()){
					gs_pt = new GlycanSequence();				
				} else if(flag == "SEQUENCE" && (next_flag == "MODIFICATION" || next_flag == "COMBINATION")) {
					// Add branch.
								
				} 
				flag = next_flag;

			} else {
				if(flag == "SEQUENCE") { // Dealing with branches.				
					//bc_pt.release();

					// Split the sequence into modules.
					sregex res = +_s;
					sregex_token_iterator cur(ROW.begin(), ROW.end(), res, -1), end;
					//Branch bc(*begin);

					size_t mono_id = 0;
					bc_pt = new Branch(boost::lexical_cast<size_t>(*cur));

					for(++cur; cur != end; cur++)
					{	
						std::string str = *cur;
						res = bos >> +_w >> eos;
						if(regex_match(str, what, res)) {
							Monosaccharide ms = fun_table.getFunctionalGroupBySymbol(str);
							bc_pt->addUnit(ms);				
						} else {
							// If linkage.
							res = bos >> (s1= +_w) >> (s2= +_d) >> '-' >> (s3= +_d) >> eos;
							if(regex_match(str, what, res)){
								Linkage lk(mono_id, boost::lexical_cast<size_t>(what[2]), boost::lexical_cast<size_t>(what[3]), boost::lexical_cast<std::string>(what[1]));
								bc_pt->addLinkage(lk);
								mono_id++;
							}
						}
					}
					gs_pt->addBranch(*bc_pt);
					delete bc_pt;	

				} else if(flag == "MODIFICATION"){
					sregex res = bos >> (s1= +_d) >> '\t' >> (s2= +_d) >> '\t' >> (s3= +_d) >> (s4= _w) >> (s5= -*~_ln); 
					if(regex_match(ROW, what, res)){

						// Modification type.
						//std::string mod_string = what[5];
						

						size_t bc_id = boost::lexical_cast<size_t>(what[1]);
						size_t mono_id = boost::lexical_cast<size_t>(what[2]);
						size_t site_id = boost::lexical_cast<size_t>(what[3]);
						std::string atom = what[4];
						std::string mod_string = what[5];

						ModificationPosition mod_pos(bc_id, mono_id);

						Modifier modifier;
						if(mod_string == "2-AB") {
							modifier.modifyGlycanSequence(*gs_pt, mod_string, mod_pos);
						} else {
							modifier.modifyGlycanSequence(*gs_pt, mod_string, mod_pos, site_id);
						}

					}

				} else if(flag == "COMBINATION"){
					sregex res = bos >> (s1= +_w) >> '\t' >> (s2= +_d);
					if(regex_match(ROW, what, res)){
						std::string mod = what[1];
						size_t num = boost::lexical_cast<size_t>(what[2]);
						mod_list.insert(make_pair(mod, num));						
					}

				}

			} 

		}
	}
	seqfile.close();
}

//vector<Peak> readPeakList()
//{
//	ifstream infile ("../../test/Hex7_6.txt");
//	string line;
//	vector<Peak> peaks;
//	if(infile.is_open()) {
//		while(getline(infile, line))
//		{
//			istringstream is;
//			is.str(line);
//			double k1; float k2;
//			is >> k1 >> k2;
//			Peak pk(k1, k2);
//			peaks.push_back(pk);
//		}
//	}
//	infile.close();
//	return peaks;
//}

MonoList readMonoList()
{
	ifstream infile ("../data/Boons_11_EDD_2-_mono.txt");
	string line;
	MonoList mono;
	if(infile.is_open()) {
		// Remove title.
		getline(infile, line);
		while(getline(infile, line))
		{
			istringstream is;
			is.str(line);
			double k1; int k2;
			is >> k1 >> k2;
			mono.insert(std::make_pair(k1, k2));
		}
	}
	infile.close();
	return mono;
}

// Find the closest TreeNode and satellite with give mass.
// Binsearch algorithm.
// The algorithm should be tailored to allow massive mass matching.
// The current position of the mass on the sr_table should be recorded.
SatelliteReverse fitLibrary(SatelliteReverse& sr_table, double ms)
{
	// Lower bound, upper bound and mid position.
	size_t lb = 0, ub = sr_table.size()-1, mid;

	double ppm1 = -3e-6;
	double ppm2 = 3e-6;

	SatelliteReverse match_results;

	for(; lb <= ub; )
	{
		//cout << lb << " " << ub << endl;
		mid = (lb+ub)/2;
		double distance = (ms - sr_table.at(mid).first)/ms;
		if(distance > ppm2) {
			lb = mid + 1;
		} else if(distance < ppm1) {
			if(mid == 0)
				break;
			ub = mid - 1;
		} else { // Mass matched!
			match_results.push_back(sr_table.at(mid));
			// Update the lower bound for next matching.
			//idx = mid+1;
			// Expand.
			// 1. To the right.
			for(size_t cur = mid+1; cur <= ub; cur++)
			{
				distance = (ms - sr_table.at(cur).first)/ms;
				if(distance <= ppm2 && distance >= ppm1) {
					match_results.push_back(sr_table.at(cur));
					// Update the lower bound for next matching.
					//idx = cur + 1;
				}
				else
					break;
			}
			// 2. To the left.
			for(size_t cur = mid-1; cur >= lb; cur--)
			{
				distance = (ms - sr_table.at(cur).first)/ms;
				if(distance <= ppm2 && distance >= ppm1)
					match_results.push_back(sr_table.at(cur));
				else
					break;
			}
			break;
		}
	}
	return match_results;
}


void matchPeakList(/*vector<Peak>& peaklist*/ MonoList& monolist, LibraryTree& libtree, std::string& filename)
{
	// Sort tree nodes and peaks in an ascending order.
	// sort(tree_vec.begin(), tree_vec.end(), TreeNode::massSmaller);
	size_t chain_length = libtree.getGlycanSequence().getBranchByID(0).getUnitNum();
	SatelliteReverse& sr_table = libtree.getReverseTable();
	//sort(peaklist.begin(), peaklist.end(), Peak::mzSmaller);
	sort(sr_table.begin(), sr_table.end(), &sateSmaller);

	//size_t charge_upper = 2;
	//double hydrogen = Composition("H").getMass();
	double electron_mass = 0.00055;
	ofstream myfile(filename.c_str());
	if(myfile.is_open()) {
		myfile << "M/Z\tZ\tInt\tM0/Z\tPPM\tAssign\n";
		
		MonoList::iterator mono_iter = monolist.begin();

		size_t mem_cur = 0;

		for(; mono_iter != monolist.end() && mem_cur < sr_table.size(); mono_iter++)
		{
			// Charge and mz.
			int k = mono_iter->second; double mz = mono_iter->first;
			double hydrogen = Composition("H").getMass();
			double mass = k * mz + k * Composition("H").getMass() - k * electron_mass;

			SatelliteReverse results = fitLibrary(sr_table, mass);

			if(results.size() == 0) // No results.
				continue;
			else {

				SatelliteReverse::iterator iter = results.begin();
				for(; iter != results.end(); iter++)
				{
					// print the results.
					//printf("%.5f", peaklist.at(m).getMZ());
					myfile << std::setprecision(7) << mz;
					myfile << "\t" << k << "\t" << "0.0\t";
					// print Cleavage type.
					TreeNodePtr tn_ptr = iter->second.first;
					const SatelliteMass& sm = tn_ptr->satellites.at(iter->second.second);
					// Theoretical m/z.
					double theo_mz = (sm.second+k*electron_mass-k*Composition("H").getMass())/k;
					//printf("%.5f", mass_charge);
					// PPM.
					myfile << std::setprecision(7) << theo_mz;
					double exp_ppm = (mz - theo_mz)/theo_mz * 1e6;
					myfile << "\t" << exp_ppm << "\t";
					const CleavageCollection& cc = tn_ptr->fg.getCleavages();
					for(CleavageCollection::const_iterator const_iter = cc.begin(); const_iter != cc.end(); const_iter++)
					{
						for(std::vector<FragmentPosition>::const_iterator const_iter2 = (*const_iter).second.begin(); const_iter2 != (*const_iter).second.end(); const_iter2++) {
							myfile << const_iter->first;
							if(const_iter->first == "A" || const_iter->first == "B" || const_iter->first == "C"){
								myfile << (*const_iter2).mono_id + 1;
							} else if(const_iter->first == "X" || const_iter->first == "Y" || const_iter->first == "Z"){
								myfile << chain_length-(*const_iter2).mono_id-1;
							}
							if(const_iter->first == "X" || const_iter->first == "A")
								myfile << ":" << const_iter2->xring_first << "-" << const_iter2->xring_second << " ";
							else 
								myfile << " ";
						}
					}
					myfile << "\t";
					//const SatelliteMass& sm = tn_ref.getSatellites().at(iter->second.second);
					Satellite::const_iterator const_iter = sm.first.begin();
					for(; const_iter != sm.first.end(); const_iter++)
					{
						// The combination of mass shift.
						if(const_iter->second != 1)
							myfile << const_iter->second;
						myfile << const_iter->first << " ";
					}
					myfile << endl;
					// Get mass.
				}
			}


		}
			//for(size_t m = 0; m < peaklist.size() && mem_cur < sr_table.size(); m++)
			//{
			//
			//	double hydrogen = Composition("H").getMass();
			//	double mass = k * peaklist.at(m).getMZ() + k * Composition("H").getMass() - k * electron_mass;

			//	SatelliteReverse results = fitLibrary(sr_table, mass, mem_cur);
			//	
			//	if(results.size() == 0) // No results.
			//		continue;
			//	else {

			//		SatelliteReverse::iterator iter = results.begin();
			//		for(; iter != results.end(); iter++)
			//		{
			//			// print the results.
			//			//printf("%.5f", peaklist.at(m).getMZ());
			//			myfile << std::setprecision(7) << peaklist.at(m).getMZ();
			//			myfile << "\t" << k << "\t" << peaklist.at(m).getIntensity() << "\t";
			//			// print Cleavage type.
			//			TreeNodePtr tn_ptr = iter->second.first;
			//			const SatelliteMass& sm = tn_ptr->satellites.at(iter->second.second);
			//			// Theoretical m/z.
			//			double mass_charge = (sm.second+k*electron_mass-k*Composition("H").getMass())/k;
			//			//printf("%.5f", mass_charge);
			//			// PPM.
			//			myfile << std::setprecision(7) << mass_charge;
			//			double exp_ppm = (peaklist.at(m).getMZ() - mass_charge)/mass_charge * 1e6;
			//			myfile << "\t" << exp_ppm << "\t";
			//			const CleavageCollection& cc = tn_ptr->fg.getCleavages();
			//			for(CleavageCollection::const_iterator const_iter = cc.begin(); const_iter != cc.end(); const_iter++)
			//			{
			//				for(std::vector<FragmentPosition>::const_iterator const_iter2 = (*const_iter).second.begin(); const_iter2 != (*const_iter).second.end(); const_iter2++) {
			//					myfile << const_iter->first;
			//					if(const_iter->first == "A" || const_iter->first == "B" || const_iter->first == "C"){
			//						myfile << (*const_iter2).mono_id + 1;
			//					} else if(const_iter->first == "X" || const_iter->first == "Y" || const_iter->first == "Z"){
			//						myfile << chain_length-(*const_iter2).mono_id-1;
			//					}
			//					if(const_iter->first == "X" || const_iter->first == "A")
			//						myfile << ":" << const_iter2->xring_first << "-" << const_iter2->xring_second << " ";
			//					else 
			//						myfile << " ";
			//				}
			//			}
			//			myfile << "\t";
			//			//const SatelliteMass& sm = tn_ref.getSatellites().at(iter->second.second);
			//			Satellite::const_iterator const_iter = sm.first.begin();
			//			for(; const_iter != sm.first.end(); const_iter++)
			//			{
			//				// The combination of mass shift.
			//				if(const_iter->second != 1)
			//					myfile << const_iter->second;
			//				myfile << const_iter->first << " ";
			//			}
			//			myfile << endl;
			//			// Get mass.
			//		}
			//	}


			//}


	} else {
		cout << "Unable to open file: " << filename << endl;
	}

}


int main(int argc, char* argv[])
{
	Param& param = Param::Instance();
	param.setParameter<std::string>("B_type_cleavage_shift", "O");
	param.setParameter<std::string>("Y_type_cleavage_shfit", "O");

	FunctionalGroupTable& fgt = FunctionalGroupTable::Instance();

	// Load Fragmentation parameters.
	FragmentationTable& ft = FragmentationTable::Instance();
	//ft.load();

	std::vector<GlycanSequence> gag_list;
	// The name of the modification and the number.
	std::map<std::string, size_t> mod_list;
	readSequence(gag_list, mod_list);

	//vector<Peak> peaklist = readPeakList();
	MonoList monolist = readMonoList();

	for(std::vector<GlycanSequence>::iterator iter = gag_list.begin(); iter != gag_list.end(); iter++)
	{
		// For each sequence. Do the cleavage.
		//chain_length = (*iter).getBranches().front().getGlycanChainUnits().size();
		//Cleavages(*iter);
		LibraryForest libforest(*iter);
		// The library forest is too huge to be stored in memory. It's better to store only the modification pattern, but leave the rest of the calculation to the last step.
		libforest.generateTreeByComposition(mod_list["Sulfation"], mod_list["Acetation"]);
		
		for(std::vector<GlobalModPattern>::iterator iter1 = libforest.getGlobalModification().begin(); iter1 != libforest.getGlobalModification().end(); iter1++)
		{
			// update sequence and get the tree.
			LibraryTree tree = libforest.applyModification(*iter1);
			//cout << "Sequence: " << tree.getGAGSequenceCodeString() << endl;
			//matchPeakList(peaklist, *iter1);
			std::string filename("../../test/data/");
			filename.append(tree.getGAGSequenceCodeString());
			filename.append(".txt");
			//matchPeakList(peaklist, tree, filename);
			matchPeakList(monolist, tree, filename);
			cout << "Sequence: " << tree.getGAGSequenceCodeString() << " processed!" << endl;
		}	
	}
	
}