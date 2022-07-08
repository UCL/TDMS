typedef struct link *Link;
typedef struct queue *Queue;

Queue initQueue();
void ql_initLink(Link l, char *text, Link next);
void ql_freeLink(Link l);
void ql_init(Queue q);
void ql_push(Queue q,char *message);
void ql_pop(Queue q, char *messptr, int mess_size);
void ql_peek(Queue q, char *messptr, int mess_size);
int ql_getNumElements(Queue q);
void ql_deleteQueue(Queue q);
