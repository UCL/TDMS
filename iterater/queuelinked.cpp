#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "queuelinked.h"
/*Implements a simple linked FIFO queue
  
  Peter Munro, TRL, Nov 2001

  Has head and tail links which have no message. When a message is
  added to the queue it is placed before the tail (ie, last before 
  tail) and when one is removed it is taken from just after the
  head:    
               Add new link here (push)
                   /
                  /
  +-+   +-+   +-+ | +-+
  |-|-->|-|-->|-|-->|-|
  +-+ | +-+   +-+   +-+
   |   \             |
  head  \           tail
         \
Remove link from here (pop)

  Modifications:

  1) 24/1/2002, Peter Munro, TRL, Added ql_deleteQueue and modified
  ql_freeLink so that the actual link structure is freed, not just the
  text.  

*/ 

struct link{ char *text; Link next; };

struct queue{
  Link head,tail;
  int numelements;
};

Queue initQueue(){
  return (Queue)malloc(sizeof(struct queue));
}

/*Initialises a link in the queue*/

void ql_initLink(Link l, char *text, Link next){
  l->text = (char *)malloc((strlen(text)+1)*sizeof(char));
  strcpy(l->text,text);
  l->next = next;
}

/*Free a link when it is no longer required*/

void ql_freeLink(Link l){
  free(l->text);
  free(l);
}

/*Initialise the queue*/

void ql_init(Queue q){
  q->numelements = 0;
  q->tail = (Link)malloc(sizeof(struct link));
  q->head = (Link)malloc(sizeof(struct link));
  ql_initLink(q->tail,(char *)"",NULL);
  ql_initLink(q->head,(char *)"",q->tail);

}

/*Add an element to the tail of the queue*/

void ql_push(Queue q,char *message){
  Link temp;
  q->numelements++;  
  temp = (Link)malloc(sizeof(struct link));
  ql_initLink(temp,(char *)"",NULL);
  q->tail->next = temp;
  free(q->tail->text);
  q->tail->text = (char *)malloc((strlen(message)+1)*sizeof(char));
  strcpy(q->tail->text,message);
  q->tail=temp;
}

/*Remove the first element of the queue*/

void ql_pop(Queue q, char *messptr, int mess_size){
  Link temp;
  if(q->numelements){
    if(mess_size > (int)strlen(q->head->next->text))
      strcpy(messptr,q->head->next->text);
    temp = q->head->next;
    q->head->next = q->head->next->next;
    ql_freeLink(temp);
    q->numelements--;
  }
  else
    strcpy(messptr,"");
}

/*Look at the first element of the queue*/

void ql_peek(Queue q, char *messptr, int mess_size){
  if(q->numelements){
    if(mess_size > (int)strlen(q->head->next->text))
      strcpy(messptr,q->head->next->text);
  }
  else
    strcpy(messptr,"");
}

/*Get the number of elements*/

int ql_getNumElements(Queue q){
  return q->numelements;
}

/*Delete a Queue and free all memory*/
void ql_deleteQueue(Queue q){
  Link next, current;
  current = q->head;
  next = q->head->next;
  if( next != q->tail)
    while(next != q->tail){
      ql_freeLink(current);
      current = next;
      next = next->next;
    }
  else{
    ql_freeLink(current);
  } 
  ql_freeLink(q->tail);
  free(q);
}






