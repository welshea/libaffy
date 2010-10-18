
/*
   strip_comments will remove any comments starting with c to the end
   of line.
*/

int strip_comments(char *p, char c)
{
  if (!p)
    return -1;
  for (; *p != '\0'; p++)
  {
    if (*p == c)
    {
      *p = '\0';
      break;
    }
  }
  return 0;
}
