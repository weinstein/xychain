#ifndef NOTIFICATION_H
#define NOTIFICATION_H

#include <pthread.h>

#include "mutex.h"

class Notification
{
 public:
   Notification() : fNotify(false) {
      pthread_mutex_init(&fMu, NULL);
      pthread_cond_init(&fCond, NULL);
   }

   ~Notification() {
      pthread_mutex_destroy(&fMu);
      pthread_cond_destroy(&fCond);
   }

   void Notify() {
      {
         MutexLock ml(&fMu);
         fNotify = true;
      }
      pthread_cond_broadcast(&fCond);
   }
      
   void ResetNotify() {
      {
         MutexLock ml(&fMu);
         fNotify = false;
      }
      pthread_cond_broadcast(&fCond);
   }

   bool IsNotified() const {
      MutexLock ml(&fMu);
      return fNotify;
   }

   void WaitForNotification() {
      MutexLock ml(&fMu);
      while (!fNotify) {
         ml.CondWait(&fCond);
      }
   }

 private:
   bool fNotify;
   mutable pthread_mutex_t fMu;
   mutable pthread_cond_t fCond;
};

#endif
