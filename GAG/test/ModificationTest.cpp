#include <GAGPL/CHEMISTRY/ModificationTable.h>
#include <GAGPL/CHEMISTRY/FunctionalGroupTable.h>
#include <iostream>

using namespace std;
using namespace gag;

void printOperation(const OperationSet& operations)
{
	boost::tuple<size_t, AtomOperation, std::string> operation_tuple;
	BOOST_FOREACH(operation_tuple, operations)
	{
		cout << "Position: " << operation_tuple.get<0>() << endl;
		cout << "Sign: " << operation_tuple.get<1>() << endl;
		cout << "FunctionalGroup TBC: " << operation_tuple.get<2>() << endl;
	}
}
void printPair(const PairSet& pairs)
{
	std::pair<AtomOperation, std::string> operation_pair;
	BOOST_FOREACH(operation_pair, pairs)
	{
		cout << "Sign: " << operation_pair.first << endl;
		cout << "FunctionalGroup TBC: " << operation_pair.second << endl;
	}
}

void printModifcicationInformation(Modification& mod)
{
	// Print Name.
	cout << "Name: " << mod.getModificationName() << endl;
	// Print shortcuts.
	cout << "Symbols:" << " " << mod.getSymbol() << endl;

	// Print reaction(s).
	cout << "Reactions:" << endl;
	for(std::vector<Reaction>::iterator iter = mod.getModificationReactions().begin(); iter != mod.getModificationReactions().end(); iter++)
	{
		cout << "Position: " << iter->position;
		cout << "On core: " << endl;
		printPair(iter->core_operation);

		for(std::multimap<std::string, OperationSet>::iterator sub_iter = iter->sub_fg_operation.begin(); sub_iter != iter->sub_fg_operation.end(); sub_iter++)
		{
			cout << "On functionalGroup: ";
			cout << sub_iter->first << "--------" << endl;
			printOperation(sub_iter->second);
			cout << endl;
		}

		cout << "On reactant: ";
		cout << iter->reactant_operation.first << "---------" << endl;
		printOperation(iter->reactant_operation.second);
		cout << endl;

	}

}

int main(int argc, char **argv)
{
	FunctionalGroupTable& ftable = FunctionalGroupTable::Instance();
	//ftable.load();

	FunctionalGroup glcn_2ab = ftable.getFunctionalGroupBySymbol("GlcN");

	ModificationTable& mod_table = ModificationTable::Instance();
	//mod_table.load();

	// Retrieve the Modification object information.
	Modification mod1 = mod_table.getModificationByName("2-AB labeling");
	printModifcicationInformation(mod1);

	Modification mod2 = mod_table.getModificationBySymbol("S");
	Modification mod3 = mod_table.getModificationBySymbol("PENA");

	printModifcicationInformation(mod2);

	ftable.appendModificationToFunctionalGroup(glcn_2ab, mod1);
	ftable.appendModificationToFunctionalGroup(glcn_2ab, mod2, 2);
	ftable.appendModificationToFunctionalGroup(glcn_2ab, mod2, 3);
	ftable.appendModificationToFunctionalGroup(glcn_2ab, mod2, 6);
	glcn_2ab.printFunctionalGroupInformation();

	FunctionalGroup glcn_pena = ftable.getFunctionalGroupBySymbol("GlcN");
	
	ftable.appendModificationToFunctionalGroup(glcn_pena, mod3, 1);
	glcn_pena.printFunctionalGroupInformation();

	return EXIT_SUCCESS;

}