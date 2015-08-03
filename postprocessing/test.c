
main()
{
int k, l, m;

for (k=0; k<13; k++) {
   printf(" %2d\n",k);
   if( k >= 9 ) break;
}
printf(" after: %2d\n",k);

for (l=0; l<13; ++l) {
   printf(" %2d\n",l);
   if( l >= 9 ) break;
}
printf(" after: %2d\n",l);

printf(" check: %2d  %2d\n",k,l);
}
