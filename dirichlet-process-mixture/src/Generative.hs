{-# LANGUAGE NoMonomorphismRestriction #-}

import Control.Monad
import System.Random.MWC.Probability

mixing :: Int -> Prob IO [Double]
mixing k = do
  a   <- inverseGamma 1 1
  symmetricDirichlet k a

params :: Int -> Double -> Double -> Prob IO [(Double, Double)]
params k muy vary = do
  l   <- normal muy (sqrt vary)
  r   <- gamma 1 (recip (sqrt vary))
  mus <- replicateM k (normal l (sqrt (recip r)))

  b  <- inverseGamma 1 1
  w  <- gamma 1 vary
  ss <- replicateM k (gamma b (recip w))

  return $ zip mus ss

fmm :: Int -> Double -> Double -> Prob IO [Double]
fmm k muy vary = do
  (mus, ss) <- fmap unzip (params k muy vary)
  pis       <- mixing k

  xs <- zipWithM normal mus (fmap (sqrt . recip) ss)

  return $ zipWith (*) pis xs

conditional n k muy vary = do
  (mus, ss) <- fmap unzip (params k muy vary)
  let fs = zipWithM normal mus (fmap (sqrt . recip) ss)
  pis       <- mixing k

  replicateM n $ do
    f <- fs
    return $ zipWith (*) pis f

main :: IO ()
main = do
  samples <- withSystemRandom . asGenIO $ \gen ->
    -- replicateM 5000 (sample (fmm 5 1 1) gen)
    sample (conditional 5000 5 1 1) gen
  let pretty = putStrLn . filter (`notElem` "[]") . show
  mapM_ pretty samples




