/********************************************************************
	created:	2012/06/27
	created:	27:6:2012   15:48
	filename: 	ModificationTable.h
	file path:	GAG\src\GAGPL\CHEMISTRY
	file base:	ModificationTable
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_MODIFICATIONTABLE_H
#define GAG_MODIFICATIONTABLE_H

#include <GAGPL/MISC/ConfigLoader.h>
#include <GAGPL/CHEMISTRY/Modification.h>
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

namespace gag 
{
	class ModificationTable: public ConfigLoader, private boost::noncopyable
	{
		typedef boost::shared_ptr<Modification> ModificationPtr;
	
	private:
		std::map<std::string, ModificationPtr> mod_by_name;
		std::map<std::string, ModificationPtr> mod_by_symbols;

		ModificationTable() {
			load();
		}

	public:

		static ModificationTable& Instance();

		void load(const std::string& filename = "./config/modification.xml");

		Modification getModificationBySymbol(const std::string& symbol);

		Modification getModificationByName(const std::string& name);


	};
}



#endif /* GAG_MODIFICATIONTABLE_H */