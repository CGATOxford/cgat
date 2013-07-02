/** class for calculating connected components
 */

#ifndef _connected_components_h
#define _connected_components_h

#include <vector>
#include <map>
#include <string>

template<class T>
class Components
{

  typedef int Index;

 public:
  // constructor
  Components();

  // destructor
  virtual ~Components();

  /** @brief add a link
      @param a token of first link
      @param a token of second link
      @return true, if the link joins two disconnected components
  */
  virtual bool add( const T & a, const T & b);

  // start from new
  virtual void reset();
  
  /** @brief get component id that a token belongs to.
      @param a token to look up
      @return component id or 0, if token is not found.
  */
  virtual int get( const T & a);

  /** @brief get id for a token.
      @param a token to look up
      @return id or 0, if token is not found.
  */
  virtual Index getIndex( const T & a);

  /** @brief get token for an id
      @param an id to look up
      @return a token
  */
  virtual const T & getToken( Index );

  /** @brief get component id that an id belongs to.
      @param a token to look up
      @return component id or 0, if id is not found.
  */
  virtual int getComponent( Index );

  /** @brief get number of tokens submitted
      @return number of tokens
  */
  virtual int getNumNodes();

 protected:

  typedef typename std::map< T, Index > MapToken2Vertex;
  typedef typename std::map< T, Index >::iterator MapToken2VertexIterator;
  
  std::vector< Index > mDad;

  MapToken2Vertex mMapToken2Vertex;

  std::vector< T > mMapVertex2Token;
};

typedef Components<const char *>CharComponents;
typedef Components<int>IntComponents;
typedef Components<std::string>StringComponents;

#endif
