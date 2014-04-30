#include <GAGPL/GAGLIBRARY/FullMap.h>
#include <GAGPL/IO/SequenceReader.h>

using namespace std;
using namespace gag;
using namespace gag::io;

int main(int argc, char **argv)
{
  try
  {
    // 1. Load the sequence.
    SequenceReader seq_reader;

    // 2. Create the library

    // 3. Create the map.

    // 4. Print the graph.
  }
  catch (exception& e)
  {
    cout << "Exception: " << e.what() << endl;
  }
}