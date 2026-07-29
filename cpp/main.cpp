// polybuilder / main.cpp  --  runnable C++ demo
// =============================================
// Mirrors run_example.py. Build with CMake (see CMakeLists.txt), then run:
//     ./polybuilder
#include <iostream>
#include "polymer_builder.hpp"
#include "monomers.hpp"

using namespace polybuilder;

static void show(const std::string& title, const RDKit::ROMol& m) {
  std::cout << "  " << title << " : "
            << PolymerBuilder::toSmiles(m) << "\n";
}

int main() {
  auto M = makeMonomers();
  auto I = makeInitiators();

  std::cout << "1. HOMOPOLYMERS (n = 4)\n";
  for (const char* name : {"DMAPS", "A3367", "NIPAM"}) {
    auto p = PolymerBuilder::capHydrogen(
        PolymerBuilder::homopolymer(*M[name], 4));
    show(std::string(name) + " x4", *p);
  }

  std::cout << "2. COPOLYMERS\n";
  show("DMAPS-alt-MMA",
       *PolymerBuilder::capHydrogen(
           PolymerBuilder::alternating(*M["DMAPS"], *M["MMA"], 2)));
  show("DMAPS-block-MMA",
       *PolymerBuilder::capHydrogen(
           PolymerBuilder::block(*M["DMAPS"], 2, *M["MMA"], 2)));

  std::cout << "3. INITIATOR END GROUPS\n";
  show("KPS-[A3367 x3]",
       *PolymerBuilder::addInitiator(
           PolymerBuilder::homopolymer(*M["A3367"], 3), *I["KPS"], "both"));
  show("AIBA-[DABCO x2]",
       *PolymerBuilder::addInitiator(
           PolymerBuilder::homopolymer(*M["DABCO"], 2), *I["AIBA"], "head"));

  std::cout << "4. DP FROM LOADING\n";
  int dp = PolymerBuilder::estimateDP(265.3, 270.3, /*wt%*/ 1.0);
  std::cout << "  A3367 with 1 wt% KPS -> DP ~ " << dp << "\n";
  return 0;
}
