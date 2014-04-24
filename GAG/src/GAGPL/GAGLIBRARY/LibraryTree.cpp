#include <GAGPL/GAGLIBRARY/LibraryTree.h>
#include <GAGPL/MATH/combination.h>
//#include <GAGPL/GLYCAN/MonosaccharideUnitTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <GAGPL/FRAGMENTATION/FragmentationTable.h>
//#include <boost/assign/ptr_map_inserter.hpp>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/make_shared.hpp>
#include <GAGPL/MATH/MassConversion.h>

namespace gag
{
	//TreeNode& TreeNode::operator =(const TreeNode& rhs)
	//{
	//	if(this != &rhs)
	//	{
	//		degree = rhs.degree;
	//		fg = rhs.fg;
	//		satellites = rhs.satellites;
	//	}
	//	return *this;
	//}
	//void LibraryTree::processFragment(Fragment& fg)
	//{

	//	if(fg.getMass() > 50.0) {
	//		Composition tmp_compo = fg.getComposition();
	//		// Create a tree node.
	//		//TreeNode tn(fg);
	//		TreeNodePtr tn(new TreeNode(fg));
	//		Satellite sl;
	//		processMassLoss(tn, sl, 0, tmp_compo);	

	//		fragmentation_tree.push_back(tn);

	//	}

	//}
	//void LibraryTree::addCleavage()
	//{
	//	Fragment fg(glyco_seq);
	//	// Select the branch id.
	//	for(size_t bc_id = 0; bc_id < glyco_seq.getBranchSize(); bc_id++)
	//	{
	//		for(int re_cl = RE_INTACT; re_cl < RE_N; re_cl++)
	//		{
	//			if(re_cl == RE_INTACT) // No RE type cleavage.
	//				this->addNRECleavage(fg, re_cl, bc_id, /*glyco_seq.getBranchByID(bc_id).getUnitNum(),*/ upper_cleavage);
	//			else { // One RE type cleavage.
	//				// Go over all mono_ids.
	//				for(size_t mono_id = 0; mono_id < glyco_seq.getBranchByID(bc_id).getUnitNum(); mono_id++)
	//				{
	//					if(re_cl == A){
	//						// The distance between the start and end of the cleavage should be at least 2.
	//						size_t ring_end = glyco_seq.getBranchByID(bc_id).getUnitByID(mono_id).getRingPosByRingID().second;
	//						for(size_t i = 0; i <= ring_end-2; i++)
	//						{
	//							for(size_t j = i+2; j<=ring_end; j++)
	//							{	
	//								if(i == 0 && j == ring_end) // 0-5 cleavage is not allowed.
	//									continue;
	//								if(j == i + 3)
	//									continue;
	//								// Get a copy of the fragment.
	//								Fragment fg1 = fg;
	//								FragmentPosition fp(bc_id,mono_id,i,j);

	//								//std::cout << "MONO ID: " << mono_id << " I: " << i << " J: " << j << std::endl;
	//								// Some types are not allowed for terminal A type.
	//								if(mono_id == glyco_seq.getBranchByID(bc_id).getUnitNum()-1) {
	//									if(i == 0 && j == 4) {
	//										continue;
	//									} else if(i == 1 && j == 3) {
	//										continue;
	//									}
	//								}
	//								// The internal cleavage doesn't work for A type.
	//								if(fg1.isInternalCleavage(fp) && mono_id != glyco_seq.getBranchByID(bc_id).getUnitNum()-1) {
	//									continue;
	//								}

	//								fg1.setFragmentation("A", fp);
	//								// fragment, current position. No other cleavage is allowed to be at mono_id.
	//								this->processFragment(fg1);
	//								if(upper_cleavage >= 1 && mono_id > 0) {
	//									this->addNRECleavageOnBranch(fg1, re_cl, bc_id, mono_id-1, 1);													
	//									if(upper_cleavage > 1 && upper_cleavage+1 <= glyco_seq.getBranchSize())
	//										this->addNRECleavage(fg1, re_cl, bc_id, /*mono_id,*/ upper_cleavage-1);
	//								} 
	//									
	//							}
	//						}
	//					} else {
	//						Fragment fg1 = fg;
	//						FragmentPosition fp(bc_id, mono_id,0,0);
	//						// Calculate the cleavage type based on shift from 'A'
	//						const char c = 'A' + (re_cl-1);
	//						std::string str(&c,1);
	//						fg1.setFragmentation(str, fp);
	//						this->processFragment(fg1);
	//						if(upper_cleavage >= 1 /*&& mono_id > 0*/) {
	//							this->addNRECleavageOnBranch(fg1, re_cl, bc_id, mono_id, 1);
	//							if(upper_cleavage > 1 && upper_cleavage+1 <= glyco_seq.getBranchSize())
	//								this->addNRECleavage(fg1, re_cl, bc_id, /*mono_id+1,*/ upper_cleavage-1);
	//						}
 //
	//							
	//					}
	//				}
	//			}		
	//		}

	//	}
	//	
	//}



	//void LibraryTree::addNRECleavageOnBranch(Fragment& fg, int re_cl, const size_t& bc_id, const size_t& mono_id, size_t times)
	//{
	//	
	//	for(int nre_cl = NRE_INTACT+1; nre_cl < NRE_N; nre_cl++)
	//	{
	//		size_t cur_id = mono_id;
	//		if((nre_cl == Y || nre_cl == Z) && (re_cl == B || re_cl == C)) {
	//			if(mono_id == 0)
	//				continue;
	//			else
	//				cur_id = mono_id -1 ;
	//		}
	//		for(size_t k = 0; k <= cur_id; k++)
	//		{
	//			
	//			size_t ring_end = glyco_seq.getBranchByID(bc_id).getUnitByID(k).getRingPosByRingID().second;

	//			//std::cout << "NRE_CL: " << nre_cl << " K: " << k;

	//			if(nre_cl == X) {
	//				for(size_t i = 0; i <= ring_end-2; i++)
	//				{
	//					for(size_t j = i+2; j<=ring_end; j++)
	//					{
	//						if(i == 0 && j == ring_end)
	//							continue;
	//						if(j == i + 3)
	//							continue;
	//						Fragment fg1 = fg;
	//						FragmentPosition fp(bc_id,k,i,j);
	//						//std::cout << " I: " << i << " J: " << j << std::endl;
	//						
	//						fg1.setFragmentation("X", fp);
	//						
	//						this->processFragment(fg1);
	//						
	//					}
	//				}
	//			} else {
	//				Fragment fg1 = fg;
	//				FragmentPosition fp(bc_id,k,0,0);
	//				const char c = 'X' + (nre_cl-1);
	//				std::string str(&c,1);
	//			
	//				fg1.setFragmentation(str,fp);
	//				this->processFragment(fg1);

	//			}
	//			// Print fragment.

	//		}
	//	}

	//}

	//void LibraryTree::processBranchesCleavage(Fragment& fg, BranchListPtr::iterator it1, BranchListPtr::iterator it2)
	//{
	//	if(it1 != it2) {
	//		for(int nre_cl = NRE_INTACT+1; nre_cl < NRE_N; nre_cl++)
	//		{
	//			for(size_t k = 0; k < (*it1)->getUnitNum(); k++)
	//			{
	//				size_t ring_end = (*it1)->getUnitByID(k).getRingPosByRingID().second;

	//				if(nre_cl == X) {
	//					for(size_t i = 0; i <= ring_end-2; i++)
	//					{
	//						for(size_t j = i+2; j<=ring_end; j++)
	//						{
	//							Fragment fg1 = fg;
	//							FragmentPosition fp(0,k,i,j);
	//							fg1.setFragmentation("X", fp);
	//							this->processBranchesCleavage(fg1, it1+1, it2);
	//						}
	//					}
	//				} else {
	//					Fragment fg1 = fg;
	//					FragmentPosition fp(0,k,0,0);
	//					const char c = 'X' + (nre_cl-1);
	//					std::string str(&c,1);

	//					fg1.setFragmentation(str,fp);
	//					this->processBranchesCleavage(fg, it1+1, it2);
	//				}
	//			}
	//		}
	//	} else {
	//		
	//		processFragment(fg);
	//	}
	//}
	//void LibraryTree::collectNRECleavage(Fragment& fg, BranchListPtr& branch_iters, BranchListPtr::iterator i1, std::vector<Branch>::iterator i2)
	//{
	//	
	//	for(std::vector<Branch>::iterator iter = *i1 +1; iter != i2; iter++)
	//	{
	//		*i1 = iter;
	//		// Generate the fragments based on information in branch_iters.
	//		this->processBranchesCleavage(fg, branch_iters.begin(), branch_iters.end());
	//		
	//		// Move the cursors.
	//		if(i1 != branch_iters.begin())
	//			collectNRECleavage(fg, branch_iters, i1-1, iter);
	//	}
	//}

	//void LibraryTree::addNRECleavage(Fragment& fg, int re_cl, const size_t& bc_id, size_t times)
	//{
	//	size_t ori_id = glyco_seq.getOriginBranchID(bc_id);
	//	if(ori_id == bc_id)
	//		this->addNRECleavageOnBranch(fg, re_cl, bc_id, glyco_seq.getBranchByID(bc_id).getUnitNum()-1, 1);
	//	else {
	//	// Generate given number of branch ids each time.
	//		for(size_t count = 0; count < times; count++)
	//		{
	//			// Initialization
	//			std::vector<Branch>& bc_group = glyco_seq.getBranches();
	//			BranchListPtr branch_iters;
	//			for(size_t i = ori_id; i <= bc_id; i++)
	//			{
	//				branch_iters.push_back(bc_group.begin()+i);
	//			}
	//			collectNRECleavage(fg, branch_iters, branch_iters.end()-1, bc_group.begin()+bc_id);

	//		}
	//	}
	//}

	//void LibraryTree::processMassLoss(TreeNodePtr tn, Satellite& sl, size_t cur, Composition& cp)
	//{	
	//	//MonosaccharideUnitTable& mut = MonosaccharideUnitTable::Instance();
	//	FragmentationTable& ft = FragmentationTable::Instance();

	//	MassLossWindow mlw = ft.getMassLoss();

	//	MassLoss& ml = mlw.at(cur);
	//	// Go over all number of mass loss.
	//	// For SO3, the upper limit is determined by the min of theoretical S and customized upperlimit.
	//	int max = abs(ml.upper);
	//	if(ml.loss_compo == Composition("SO3")) {
	//		std::map<std::string, int>::const_iterator iter = cp.get().find("S");
	//		if(iter != cp.get().end()) {
	//			if(max > iter->second) {
	//				max = iter->second;
	//			}
	//		} else
	//			max = 0;
	//	}

	//	//TreeNode tmp_tn(tn);
	//	for(int i = ml.lower; i <= max; i++)
	//	{
	//		// Copy temporary objects.
	//		Composition tmp_compo = cp;

	//		Satellite tmp_sl(sl);

	//		// e.g. 2H2O
	//		if(i != 0) {
	//			//cout << ml.loss_compo.getCompositionString() << ":" << i << " ";
	//			tmp_sl.insert(make_pair(ml.loss_compo.getCompositionString(), i));
	//		}

	//		if(i < 0){
	//			for(int j = 0; j < abs(i); j++) {					
	//				tmp_compo.add(ml.loss_compo);
	//			}
	//		} else if(i > 0) {
	//			// The evaluation of mass might not be correct.
	//			//if((tmp_compo.getMass() - i * ml.loss_compo.getMass()) < 0.0)
	//				//break;
	//			for(int j = 0; j < i; j++) {
	//				//std::cout << j << " X " << tmp_compo.getCompositionString() << std::endl;
	//				if(ml.loss_compo < tmp_compo)
	//					tmp_compo.deduct(ml.loss_compo);
	//				else 
	//					break;
	//			}
	//		}
	//		if(cur == mlw.size()-1) {
	//			// The end. Print out the composition and mass.
	//			//cout << tmp_compo.getCompositionString() << " ";
	//			//cout << setprecision(6) << tmp_compo.getMass() << "\n  ";

	//			tn->satellites.push_back(make_pair(tmp_sl, tmp_compo.getMass()));
	//			//sr_table.insert(tmp_compo.getMass(), make_pair(&tn, tmp_sl));

	//		} else {
	//			processMassLoss(tn, tmp_sl, cur+1,tmp_compo);
	//		} 

	//	}

	//}
 // // Basically, this is the mass reference table.
	//void LibraryTree::createReverseTable()
	//{
	//	//using namespace boost::assign;
	//	std::vector<TreeNodePtr>::iterator iter = fragmentation_tree.begin();
	//	//for(; iter != fragmentation_tree.end(); iter++)
	//	for(size_t idx = 0; idx <fragmentation_tree.size(); idx++)
	//	{
	//		for(size_t i=0; i< fragmentation_tree[idx]->satellites.size(); i++)
	//		{
	//			sr_table.insert(std::make_pair(fragmentation_tree[idx]->satellites.at(i).second, std::make_pair(fragmentation_tree[idx], i)));
	//			/*boost::ptr_map<size_t, TreeNode> node_sate;
	//			node_sate.insert(i, &(*iter));
	//			sr_table.push_back(std::make_pair(iter->satellites.at(i).second, node_sate));*/

	//		}
	//	}
	//}
	//void LibraryTree::build()
	//{
	//	// Start from the glycan sequence, sequentially add new tree nodes.
	//	this->addCleavage();
	//	this->createReverseTable();

	//}


	////std::vector<Composition> matchLibrary(double mass, double ppm);
	
	NodeItem::NodeItem(const NodeItem& node)
	{
		sate = node.sate;
		fg = node.fg;
		compo = node.compo;
	}
	NodeItem& NodeItem::operator=(const NodeItem& node)
	{
		if(this != &node)
		{
			sate = node.sate;
			fg = boost::make_shared<Fragment>(*(node.fg));
			compo = node.compo;
		}
		return *this;
	}

	bool NodeItem::appendMassShift(const std::string& loss_type, int num)
	{	
		// Update the composition information.
		Composition append_compo;
		for(int i = 1; i <= abs(num); i++)
		{
			append_compo.add(loss_type);
		}

		//Composition append_compo(loss_type);
		if(num > 0) {
			if(append_compo < compo)
				compo.deduct(append_compo);
			else
				return false;
		} else {
			compo.add(append_compo);
		}
		// Store the mass shift information.
		sate.insert(std::make_pair(loss_type, num));
		return true;
	}

	std::string NodeItem::getCompositionShift() const
	{
		std::string shift_string;
		
		for(Satellite::const_iterator iter = sate.begin(); iter != sate.end(); iter++)
		{	
			// The inforamtion of so3 and ac is recorded separately.
			if(iter->first == "O3S" || iter->first == "CH2CO")
				continue;

			if(iter->second == 0) 
				continue;
			else if(iter->second >= 1) {
				shift_string.append("-");
				if(iter->second == 1)
					shift_string.append(iter->first);
				else {
					std::ostringstream ostr;
					ostr << iter->second << iter->first;
					shift_string.append(ostr.str());
				}
			}	else if(iter->second <= -1) {
				shift_string.append("+");
				if(iter->second == -1)
					shift_string.append(iter->first);
				else {
					std::ostringstream ostr;
					ostr << abs(iter->second) << iter->first;
					shift_string.append(ostr.str());
				}

			}
		}
		return shift_string;
	}

	void NodeItem::printNodeItem(const std::string format, const size_t clv_num) const
	{
		// The function is used to control cleavage numbers.
		if(this->getCleavageNum() > clv_num)
			return;
		
		if(format == "regular") {
			std::cout << "Cleavage type: " << this->getCleavageType() << std::endl;;
			std::cout << "Shift: " << this->getCompositionShift() << std::endl;
			std::cout << "Composition: " << this->getCompositionString() << std::endl;
			std::cout << "Ac: " << this->getModificationNum("Ac") << std::endl;
			std::cout << "SO3: " << this->getModificationNum("SO3") << std::endl;
			std::cout << "Mass: " << this->getMass() << std::endl;
		} else if(format == "compact") {
			std::cout << this->getCleavageType() << "\t" << this->getCompositionShift() << "\t" << this->getCompositionString() << "\t" << this->getModificationNum("Ac") << "\t" << this->getModificationNum("SO3") << "\t" << this->getMass() << std::endl;
		}

	}

	int NodeItem::getModificationNum( const std::string& compo_string ) const
	{
		std::string internal_string;
		if(compo_string == "Ac")
			internal_string = "C2H2O";
		else if(compo_string == "SO3")
			internal_string = "O3S";
		Satellite::const_iterator iter = sate.find(internal_string);
		return iter != sate.end() ? abs(iter->second) : 0;
	}

	void FragmentTree::build()
	{
		this->addCleavage();
	}

	void FragmentTree::addCleavage()
	{
		// Get the number of mono units.
		size_t mono_num = glycan_seq->getBranchByID(0).getUnitNum();

		for(int re_clv = RE_INTACT; re_clv <RE_N; re_clv++)
		{
			
			if(re_clv == RE_INTACT) {
				//FragmentPtr fg(new Fragment(glycan_seq));
				//std::cout << "INTACT" << std::endl;
				CleavageCollection cc;
				this->addRestCleavage(cc, mono_num);
			} else {
				for(size_t i = 0; i<mono_num; i++) {
					if(re_clv == A) {
						for(size_t x = 0; x<=3; x++)
							for(size_t y = x+2; y <= 5; y++)
							{
								if((y - x == 3) || (x == 0 && (y == 4 || y == 5)) || (x == 1 && y == 3))
									continue;
								//std::cout << "A" << i << ":" << x << "-" << y << std::endl;
								//FragmentPtr fg(new Fragment(glycan_seq));
								//fg->setFragmentation("A", i, "", x, y);
								CleavageCollection cc;
								cc.insert(std::make_pair("A", FragmentPosition(0, i, x, y)));
								this->addRestCleavage(cc, i);
							}
					} else {
						//FragmentPtr fg(new Fragment(glycan_seq));

            // B and C type of cleavage at the RE is a problem.
            // Ignore it at this moment.
            if(i == mono_num-1) continue;

						std::string type;
						if(re_clv == B) 
							type = "B";
						else if(re_clv == C) 
							type = "C";
						else 
							throw std::runtime_error("Uncharacterized cleavage type");
							
						CleavageCollection cc;
						cc.insert(std::make_pair(type, FragmentPosition(0, i, 0, 0)));
						//std::cout << type << i << std::endl;
						//fg->setFragmentation(type, i, "", 0, 0);
						this->addRestCleavage(cc, i);
					}
				}
			}
			
	
		}
	}

	void FragmentTree::addRestCleavage(const CleavageCollection& cc, const size_t& cur )
	{
		// Get the number of mono units.
		size_t mono_num = glycan_seq->getBranchByID(0).getUnitNum();
		MassLossWindow mlw = frag_table.getMassLoss();
		for(int nre_clv = NRE_INTACT; nre_clv < NRE_N; nre_clv++)
		{
			// Intact sequence.
			if(nre_clv == NRE_INTACT) {
				// Do nothing.
				//std::cout << "INTACT" << std::endl;
				this->processFragment(cc, mlw);
			} else {
				for(size_t i = 0; i < cur; i++)
				{
					if(nre_clv == X) {
						for(size_t x=0; x<=3; x++)
							for(size_t y = x+2; y <=5; y++)
							{
								if((y - x == 3) || (x == 0 && (y == 4 || y == 5)) || (x == 1 && y == 3))
									continue;
								//std::cout << "X" << i << ":" << x << "-" << y << std::endl;
								//FragmentPtr fg1(new Fragment(*fg));
								//fg1->setFragmentation("X", i, "", x,y);
								CleavageCollection cc1(cc);
								cc1.insert(std::make_pair("X", FragmentPosition(0, i, x, y)));
								this->processFragment(cc1, mlw);
							}
					} else {
						//FragmentPtr fg1(new Fragment(*fg));

						std::string type;
						if(nre_clv == Y) 
							type = "Y";
						else if(nre_clv == Z) 
							type = "Z";
						else 
							throw std::runtime_error("Uncharacterized cleavage type");
						
						//std::cout << type << i << std::endl;
						//fg1->setFragmentation(type, i, "", 0, 0);
						CleavageCollection cc1(cc);
						cc1.insert(std::make_pair(type, FragmentPosition(0, i, 0, 0)));
						this->processFragment(cc1, mlw);
					}
				}
			}
		}
	}

	void FragmentTree::processFragment(const CleavageCollection& cc, MassLossWindow mlw)
	{	
		FragmentPtr fg(new Fragment(glycan_seq));
		fg->updateCleavage(cc);

		//std::string final_type = fg->getCleavageType();
		
		// Modify the MassLossWindow based on the possibility of modification sites.
		std::set<std::string> mod_types = glycan_seq->getModificationTypes();

		for(std::set<std::string>::iterator iter = mod_types.begin(); iter != mod_types.end(); iter++)
		{
			int num = (int)fg->getModificationSiteNum(*iter, 1);
			Composition temp_compo;
			if(*iter == "SO3") {
				temp_compo.add("SO3");
				int upper_num = glycan_seq->getModificationConstraint("SO3");
				if(num > upper_num)
					num = upper_num;
			} else if(*iter == "Ac") {
				temp_compo.add("CH2CO");
				int upper_num = glycan_seq->getModificationConstraint("Ac");
				if(num > upper_num)
					num = upper_num;
			}
			if(num != 0) mlw.push_back(MassLoss(temp_compo, 0, -1 * num));

		}
		
		this->processFragment(fg, mlw, 0, Satellite());
	}

	void FragmentTree::processFragment( FragmentPtr fg, MassLossWindow& mlw, size_t degree, Satellite sate )
	{
		if(degree == mlw.size()) {  // Out of boundary.
			NodeItem node(fg);

			int max_num = (int)fg->getModificationSiteNum("SO3", 1);
			Satellite::iterator ac_iter = sate.find("CH2CO");
			int ac_num = ac_iter == sate.end() ? 0 : -1 * ac_iter->second;
			Satellite::iterator s_iter = sate.find("SO3");
			int s_num = s_iter == sate.end() ? 0 : -1 * s_iter->second;
			if(s_num + ac_num > max_num)
				s_iter->second = max_num - ac_num;

			// Update the composition.
			for(Satellite::iterator iter = sate.begin(); iter != sate.end(); iter++) {
				//std::cout << iter->first << " " << iter->second << " ";
				if(!node.appendMassShift(iter->first, iter->second)) {
					return;
				}
			}
			//std::cout << std::endl;

			node_pool.insert(node);
			return;
		}

		MassLoss& ml = mlw.at(degree);

		for(int i = ml.lower; i<=ml.upper; i++)
		{
			// Copy the composition.
			Satellite temp_sate(sate);

			if(i != 0) {
				temp_sate.insert(std::make_pair(ml.loss_compo.getCompositionString(), i));
				this->processFragment(fg, mlw, degree+1, temp_sate);
			} else {
				this->processFragment(fg, mlw, degree+1, temp_sate);
			}
		}

	}

	void FragmentTree::getNodeItemsByMass(MonoPeakPtr pk, std::multimap<MonoPeakPtr, NodeItem>& node_map)
	{
		TreeByMass& tree_mass_index = node_pool.get<theo_mass>();
		double matching_error = param.getParameter<double>("matching_error").first;

		double mass = msmath::calculateMass(pk->mz, -1 * pk->z);
		TreeByMass::iterator mass_iter = tree_mass_index.lower_bound(mass);
    TreeByMass::iterator up_bound = tree_mass_index.upper_bound(mass);
		// Shift the iterator upstream and downstream.
		
		if(mass_iter == tree_mass_index.end() && up_bound == tree_mass_index.end()) // Not found.
			return;

		//std::cout << "Current mass: " << mass << std::endl;
		TreeByMass::iterator up_iter = mass_iter;

		// Upstream.
		while(1)
		{	
			//std::cout << "Try to match " << up_iter->getMass() << std::endl;
			if((up_iter->getMass() - mass)/mass >= 1e-6 * matching_error)
				break;

			if(abs(up_iter->getMass()-mass)/mass < 1e-6 * matching_error){
				//std::cout << "A hit for " << mass << std::endl;
				node_map.insert(std::make_pair(pk, *up_iter));	
				//if(status == false) status = true;
			} 
			up_iter++;
			if(up_iter == tree_mass_index.end())
				break;
		}
		TreeByMass::iterator down_iter = mass_iter;
		while(1)
		{
			//std::cout << "Try to match " << down_iter->getMass() << std::endl;
			if((mass - down_iter->getMass())/down_iter->getMass() >= 1e-6 * matching_error) {
				break;
			}

			if(abs(down_iter->getMass()-mass)/down_iter->getMass() < 1e-6 * matching_error){
				//std::cout << "A hit for " << mass << std::endl;
				node_map.insert(std::make_pair(pk, *down_iter));
				//if(status == false) status = true;
			} 
			if(down_iter == tree_mass_index.begin())
				break;
			down_iter--;
		}

		return;
		
	}

	void FragmentTree::printLibrary()
	{
		TreeByMass& tree_mass_index = node_pool.get<theo_mass>();
		for(TreeByMass::const_iterator iter = tree_mass_index.begin(); iter!= tree_mass_index.end(); iter++)
		{
			(*iter).printNodeItem("compact");
		}
	}

	std::multimap<MonoPeakPtr, NodeItem> FragmentTree::searchLibrary( const std::set<MonoPeakPtr>& pk_list )
	{
		std::multimap<MonoPeakPtr, NodeItem> node_list;
		// Slow version.
		for(std::set<MonoPeakPtr>::const_iterator const_iter = pk_list.begin(); const_iter != pk_list.end(); const_iter++)
		{
      double mz = (*const_iter)->mz;
			this->getNodeItemsByMass(*const_iter, node_list);
		}
		return node_list;
	}

	void FragmentTree::exportLibrary( const std::string& filename )
	{
		std::ofstream outfile(filename.c_str());

		if(outfile.is_open()) {
			TreeByMass& tree_mass_index = node_pool.get<theo_mass>();
			outfile.precision(5);
			for(TreeByMass::const_iterator iter = tree_mass_index.begin(); iter!= tree_mass_index.end(); iter++)
			{
				outfile << std::fixed << iter->getMass() << "\t" << iter->getCleavageType() << "\t" << iter->getCompositionShift() << "\t" << iter->getCompositionString() << "\t" << "\n";
			}
			outfile.close();
		} else {
			std::cout << "Unable to open file!\n" << std::endl;
		}
	}

	void loadMonoPeakList( const std::string& filename, std::set<MonoPeakPtr>& peak_list)
	{
		std::cout << "Load peak list..." << std::endl;
		std::ifstream infile(filename.c_str());

		std::string line;

		//RichList spec;
		if(infile.is_open())
		{
			// Deal with title.
			std::getline(infile, line);
/*			std::istringstream title;
			title.str(line);
			double t1, t2;
			title >> t1 >> t2;		*/	

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
				peak_list.insert(pk);
			}
		}
		infile.close();
		infile.clear();

	}

}