-- stripspaces.hs

-- Thomas Miedema
-- March 2010

-- Input: filenames
-- Action: Strip all spaces between parenthesis and brackets from files,
--         without messing up formatting.

module Main where

import Text.ParserCombinators.Parsec
import Text.ParserCombinators.Parsec.Token
import Monad
import System( getArgs, system )
import Directory( doesFileExist )
import System.IO


-------------------- MAIN ------------------------

main = do
  files <- getArgs
  mapM_ stripSpacesFromFile files

-- Read from stdin, write to stdout.
--main = interact stripSpacesFromString


-------------------- STRIP SPACES ------------------------

-- Read file, strip spaces between parenthesis and brackets, overwrite file.
stripSpacesFromFile :: String -> IO()
stripSpacesFromFile file = do
    exists <- doesFileExist file
    when exists $ do input <- readFile file
		     let output = stripSpacesFromString input
		     forceList output `seq` writeFile file output

-- Strip spaces between parenthesis and brackets from string.
stripSpacesFromString :: String -> String
stripSpacesFromString input =
  case parse (pythonCode anyCharAsString 0) "" input of
    --Todo: Left err -> err
    Right output -> output

-- Match Python code.
-- The first character is either the beginning of a parens term, a brackets 
-- term, or just any character (the starting charParser is anyCharAsString).
-- The integer `level` counts the level of nestings.
pythonCode :: Parser [Char] -> Int -> Parser [Char]
pythonCode charParser level = do
  termOrChar <- try (parens2 level) <|> try (brackets2 level) <|> charParser
  rest <- pythonCode charParser level
  return $ termOrChar ++ rest
  <|> return ""

-- Match parentheses and content in between.
parens2 = term ('(',')')

-- Match brackets and content in between.
brackets2 = term ('[',']')

-- Match cOpen, content, cClose.
-- Remove `level` spaces after a newline if a newline is found.
term :: (Char, Char) -> Int -> Parser [Char]
term (cOpen, cClose) level = do
  (open, newLevel) <- opener cOpen level
  let newCharParser = newlineOrAnyCharAsStringExcept (closer cClose) newLevel
  content <- pythonCode newCharParser newLevel
  close <- closer cClose
  return $ open:[] ++ content ++ close

-- Match character cOpen, remove at most one space and count how many.
opener :: Char -> Int -> Parser (Char, Int)
opener cOpen level = do char cOpen
			n_spaces <- option 0 (char ' ' >> return 1)
			return (cOpen, level + n_spaces)

-- Remove spaces and match character cClose.
closer :: Char -> Parser [Char]
closer cClose = do beforeClose <- closeOptions
		   char cClose
		   return $ beforeClose ++ cClose:[]
  
-- Match some edge cases.
-- Don't remove spaces before single closer after a newline.
-- Keep one space between a comma and the character cClose [1, 2, 3, ].
closeOptions = do char '\n'
		  spaces <- many $ char ' '
		  return $ '\n':spaces
	   <|> do char ','
		  skipMany $ char ' '
		  return ", "
	   <|> do skipMany $ char ' '
		  return ""

-- Match newline or any char as string except exclude.
-- Remove `level` spaces after newline.
newlineOrAnyCharAsStringExcept :: Parser [Char] -> Int -> Parser [Char]
newlineOrAnyCharAsStringExcept exclude level = do notFollowedBy' exclude
						  char '\n'
						  count level $ char ' '
						  return "\n"
					   <|> anyCharAsStringExcept exclude

-- Match any char as string except exclude.
anyCharAsStringExcept :: Parser [Char] -> Parser [Char]
anyCharAsStringExcept exclude = do notFollowedBy' exclude
				   anyCharAsString

-- Match any char, return as string.
anyCharAsString :: Parser [Char]
anyCharAsString = do c <- anyChar
		     return $ c:[]


-------------------- HELPERS ------------------------

-- Strict evaluation of list.
-- Could use System.IO.Strict instead, or Data.ByteString.
forceList [] = []
forceList (x:xs) = forceList xs `seq` (x:xs)


-- KEEP THIS.
-- http://www.mail-archive.com/haskell-cafe@haskell.org/msg10552.html
-- The problem is that notFollowedBy has type
--
-- notFollowedBy  :: Show tok => GenParser tok st tok -> GenParser tok st ()
--
-- ie, the result type of the parser you pass to notFollowedBy has to be
-- the same as the token type, in this case Char.  
--
-- A solution to this found here:
-- http://www.haskell.org/pipermail/haskell/2004-February/013632.html
notFollowedBy' :: Show a => GenParser tok st a -> GenParser tok st ()
notFollowedBy' p  = try $ join $  do  a <- try p
				      return (unexpected (show a))
				  <|>
				  return (return ())

-- Todo: parens and brackets are also defined here, but a bit different:
-- http://hackage.haskell.org/packages/archive/parsec/3.0.0/doc/html/Text-ParserCombinators-Parsec-Token.html


-------------------- TESTS ------------------------

test = mapM_ test' [-- Parenthesis.
		    ("foo( bar )baz", "foo(bar)baz"),

		    -- Brackets.
		    ("foo[ bar ]baz", "foo[bar]baz"),

		    -- Unbalanced spaces.
		    ("foo[ bar]baz", "foo[bar]baz"),
		    ("foo[bar ]baz", "foo[bar]baz"),

		    -- Only one, do nothing.
		    ("foo[ bar", ""),
		    ("foo ]bar", ""),

		    -- Not matching.
		    ("foo[ bar ] ]baz", "foo[bar] ]baz"),
		    ("foo[ [ bar ]baz", "foo[ [bar]baz"),

		    -- Empty lists.
		    ("foo[ ]", "foo[]"),
		    ("foo[  ]", "foo[]"),

		    -- Nested.
		    ("foo[ [ bar ] ]baz", "foo[[bar]]baz"),

		    -- Nested and different.
		    ("foo( [ bar ] )baz", "foo([bar])baz"),

		    -- Mix. It happens to ignore the first.
		    ("foo( [ bar ) ]baz", "foo( [bar )]baz"),
		    ("foo[ ( bar ] )baz", "foo[ (bar ])baz"),

		    -- Final character is a comma of list.
		    ("foo[ 1,2,3,  ]", "foo[1,2,3, ]"),
		    ("foo[,]", "foo[, ]"),

		    -- Remove space after newline within term.
		    ("foo[ bar\n     baz ]", "foo[bar\n    baz]"),

		    -- Remove two spaces after newline within nested term.
		    ("foo[ [ bar\n       baz ] ]", "foo[[bar\n     baz]]"),

		    -- Keep space after newline outside of term.
		    ("foo\n bar", ""),
		    ("foo[ bar\n baz", ""),

    -- From epdp/utils.py.
    -- Strip at most 1 space after opener to preserve layout.
     ("M = numpy.array( [ [    0.0, - a[2],   a[1] ],\n\
      \                   [   a[2],    0.0, - a[0] ],\n\
      \                   [ - a[1],   a[0],    0.0 ] ] )", 
      "M = numpy.array([[   0.0, - a[2],   a[1]],\n\
      \                 [  a[2],    0.0, - a[0]],\n\
      \                 [- a[1],   a[0],    0.0]])"),

    -- Remove three spaces after newline within 3x nested term.
    ("M = numpy.array( [ [ cosalpha + cosalphac * r[0] * r[0],\n\
     \                     cosalphac * r[0] * r[1] - r[2] * sinalpha,\n\
     \                     cosalphac * r[0] * r[2] + r[1] * sinalpha ],\n\
     \                   [ cosalphac * r[0] * r[1] + r[2] * sinalpha,\n\
     \                     cosalpha + cosalphac * r[1] * r[1],\n\
     \                     cosalphac * r[1] * r[2] - r[0] * sinalpha ],\n\
     \                   [ cosalphac * r[0] * r[2] - r[1] * sinalpha,\n\
     \                     cosalphac * r[1] * r[2] + r[0] * sinalpha,\n\
     \                     cosalpha + cosalphac * r[2] * r[2] ] ] )",
     "M = numpy.array([[cosalpha + cosalphac * r[0] * r[0],\n\
     \                  cosalphac * r[0] * r[1] - r[2] * sinalpha,\n\
     \                  cosalphac * r[0] * r[2] + r[1] * sinalpha],\n\
     \                 [cosalphac * r[0] * r[1] + r[2] * sinalpha,\n\
     \                  cosalpha + cosalphac * r[1] * r[1],\n\
     \                  cosalphac * r[1] * r[2] - r[0] * sinalpha],\n\
     \                 [cosalphac * r[0] * r[2] - r[1] * sinalpha,\n\
     \                  cosalphac * r[1] * r[2] + r[0] * sinalpha,\n\
     \                  cosalpha + cosalphac * r[2] * r[2]]])"),

     -- Edge case from epdp/myrandom.py.
     -- Don't remove spaces before single closer after a newline.
    ("__all__ = (\n\
     \    'shuffle',\n\
     \    'uniform',\n\
     \    'normal',\n\
     \    'seed',\n\
     \    'get_raw'\n\
     \    )", "")]

  where
    test' :: (String, String) -> IO()
    test' (a, b) = let a' = stripSpacesFromString a in 
		     if (a' == a) || (a' == b)
		       then print True
		       else do putStrLn "\nInput:"
			       putStrLn a
			       putStrLn "\nExpected output:"
			       putStrLn b
			       putStrLn "\nOutput:"
			       putStrLn a'
			       putStrLn "\n"


