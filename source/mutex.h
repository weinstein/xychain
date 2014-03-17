#ifndef MUTEX_H
#define MUTEX_H

#include <pthread.h>

class MutexLock {
 public:
   void Lock() {
      if (fMu) {
         pthread_mutex_lock(fMu);
      }
   }

   void Unlock() {
      if (fMu) {
         pthread_mutex_unlock(fMu);
      }
   }

   MutexLock(pthread_mutex_t* mu) : fMu(mu) {
      Lock();
   }

   ~MutexLock() {
      Unlock();
   }

   pthread_mutex_t operator*() {
      return *fMu;
   }

   pthread_mutex_t* operator->() {
      return fMu;
   }

   pthread_mutex_t* get() {
      return fMu;
   }

   pthread_mutex_t* Release() {
      Unlock();
      pthread_mutex_t* ret = fMu;
      fMu = NULL;
      return ret;
   }

   void Reset(pthread_mutex_t* new_mu) {
      Unlock();
      fMu = new_mu;
      Lock();
   }

   void CondWait(pthread_cond_t* cond) {
      pthread_cond_wait(cond, fMu);
   }

 private:
   pthread_mutex_t* fMu;
};

#endif
