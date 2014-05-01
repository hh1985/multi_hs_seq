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