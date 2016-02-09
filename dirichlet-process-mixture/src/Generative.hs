
-- | A generative (prior predictive) finite mixture model.

module Generative where

import Control.Monad (replicateM, zipWithM)
import System.Random.MWC.Probability

-- | Mixing probabilities for the model.
mixing :: Int -> Prob IO [Double]
mixing k = do
  a <- inverseGamma 1 1
  symmetricDirichlet k a

-- | Mean and precision parameters for clusters in the model.
locationScale :: Int -> Double -> Double -> Prob IO [(Double, Double)]
locationScale k muy vary = do
  l   <- normal muy (sqrt vary)
  r   <- gamma 1 (recip (sqrt vary))
  mus <- replicateM k (normal l (sqrt (recip r)))
  b   <- inverseGamma 1 1
  w   <- gamma 1 vary
  ss  <- replicateM k (gamma b (recip w))
  return $ zip mus ss

-- | Finite Gaussian mixture model.
fmm :: Int -> Double -> Double -> Prob IO [Double]
fmm k muy vary = do
  (mus, ss) <- fmap unzip (locationScale k muy vary)
  pis <- mixing k
  xs  <- zipWithM normal mus (fmap (sqrt . recip) ss)
  return $ zipWith (*) pis xs

-- | Conditional mixture model.  Sample a set of probabilities from the mixing
--   model and then use those to generate observations in each cluster.
conditional :: Int -> Int -> Double -> Double -> Prob IO [[Double]]
conditional n k muy vary = do
  (mus, ss) <- fmap unzip (locationScale k muy vary)
  pis <- mixing k
  replicateM n $ do
    f <- zipWithM normal mus (fmap (sqrt . recip) ss)
    return $ zipWith (*) pis f

-- | Sample 5000 times from a conditional mixture model.
main :: IO ()
main = do
  samples <- withSystemRandom . asGenIO $ sample (conditional 5000 5 1 1)
  let pretty = putStrLn . filter (`notElem` "[]") . show
  mapM_ pretty samples

