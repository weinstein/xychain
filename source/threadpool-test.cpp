#include "threadpool.h"

#include <assert.h>
#include <iostream>
#include <unistd.h>
using namespace std;

#include "notification.h"
#include "testing.h"

int main(int argc, char ** argv) {
   return 0;
}

void do_work(void* arg) {
   bool* completed = (bool*)arg;

   sleep(1);

   *completed = true;
}

DEFINE_TEST(testExecuteAndWait) {
   ThreadPool workpool(4);

   bool completed[4] = {false, false, false, false};
   for (int i = 0; i < 4; i++) {
      workpool.AddJob(do_work, (void*)(&completed[i]));
   }
   workpool.Wait();

   for (int i = 0; i < 4; i++) {
      assert(completed[i] == true);
   }
}

DEFINE_TEST(testNotification) {
   ThreadPool workpool(2);

   bool completed[4] = {false, false, false, false};
   Notification notifiers[4];
   for (int i = 0; i < 4; i++) {
      workpool.AddJobWithNotification(do_work, (void*)(&completed[i]), &notifiers[i]);
   }
   for (int i = 0; i < 4; i++) {
      notifiers[i].WaitForNotification();
      assert(notifiers[i].IsNotified());
      assert(completed[i]);
   }
   assert(workpool.GetPendingJobs() == 0);
}
