/********************************************************************
	created:	2012/06/24
	created:	24:6:2012   16:35
	filename: 	SequenceReader.h
	file path:	GAG\src\GAGPL\IO
	file base:	SequenceReader
	file ext:	h
	author:		Han Hu
	
	purpose:	
*********************************************************************/

#ifndef GAG_SEQUENCEREADER_H
#define GAG_SEQUENCEREADER_H

#include <boost/ptr_container/ptr_vector.hpp>
#include <GAGPL/GLYCAN/GlycanSequence.h>
#include <GAGPL/GAGLIBRARY/SequenceSpace.h>
#include <GAGPL/SPECTRUM/PeakList.h>
#include <GAGPL/SPECTRUM/MonoPeak.h>
#include <GAGPL/MISC/Param.h>

using namespace param;

namespace gag 
{
	namespace io
	{
		class SequenceReader 
		{
		public:
			// Constructor.
			SequenceReader()
				: param(Param::Instance())
			{}

			void readSequenceStructure(const std::string& filename, GlycanSequencePtr seq);
      // TBD: New function for read sequence.
      GlycanSequencePtr readSequenceStructure(const std::string& filename, const std::string path);

			void readSpectrum(const std::string& filename, RichList& specrum);

			void readMonoPeakList(const std::string& filename, std::set<MonoPeakPtr>& mono_list);

			//void writePredictionResults(const std::string& filename, std::string title, std::string seq_code, std::string m_code, const std::vector<double>& seq_scores);

		private:
			Param& param;
		};

	}
	
}

#endif /* GAG_SEQUENCEREADER_H */ 
