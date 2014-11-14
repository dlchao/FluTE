/* class R0model
 * Subclass of model for determining R0
 * Only the index case is infectious
 * Secondary cases are not infected; instead, they are no longer susceptible
 */

using namespace std;

class R0Model : public EpiModel {
 public:
  R0Model(EpiModelParameters &params);
  virtual void prerun(void);
  virtual void run(void);
 protected:
  virtual void infect(Person& p);
  virtual void summary(void);
  int nNumInfected;  // total number of infected people (including index)
  unsigned int nIndexCaseID;  // ID of the first person infected in the simulation
};
