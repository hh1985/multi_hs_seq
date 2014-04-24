/*
 * =====================================================================================
 *
 *       Filename:  FragmentationTable.h
 *
 *    Description:  Functional group table loaded from xml drive.
 *
 *        Version:  1.0
 *        Created:  04/29/2012  9:51:00 PM
 *       Revision:  none
 *       Compiler:  msvc
 *
 *         Author:  Han Hu (HH), hh.earlydays@gmail.com
 *   Organization:  Boston University
 *
 * =====================================================================================
 */


#ifndef	GAG_FRAGMENTATIONTABLE_H
#define GAG_FRAGMENTATIONTABLE_H

#include <GAGPL/MISC/ConfigLoader.h>
#include <GAGPL/CHEMISTRY/Composition.h>
#include <GAGPL/FRAGMENTATION/FragmentationParams.h>
#include <boost/noncopyable.hpp>

namespace gag
{
	typedef std::pair<int, Composition> CompositionSigned;
	typedef std::vector<CompositionSigned> CompositionShift;
	struct MassLoss
	{
		Composition loss_compo;
		int upper;
		int lower;

		MassLoss(const Composition& compo, int u, int l)
			: loss_compo(compo), upper(u), lower(l)
		{
		}
		MassLoss() {}
	};
	typedef std::vector<MassLoss> MassLossWindow;

	class FragmentationTable: public ConfigLoader, private boost::noncopyable
	{
		private:
			// string here is type.
			std::map<std::string, FragmentationParams> fragmentation_params;
			// string here is the signed composition
			std::map<std::string, std::pair<int, int> > mass_loss; 
		protected:
			FragmentationTable() {
				load();
			}

		public:
			
			static FragmentationTable& Instance();

			void load(const std::string& filename = "./config/fragmentation.xml");

			FragmentationParams getFragmentationParams(const std::string& cleavage_type) const;
			CompositionSigned getCompositionSigned(const std::string& str) const;
			CompositionShift getCleavageShift(const std::string& cleavage_type) const;
			CompositionShift getCleavageShift(const std::string& cleavage_type, const std::string& dis_type) const;			 
			MassLossWindow getMassLoss() const;
	};
}



#endif   /* ----- #ifndef GAG_FUNCTIONALGROUPTABLE_H ----- */
